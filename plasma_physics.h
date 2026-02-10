#ifndef PLASMA_PHYSICS_H
#define PLASMA_PHYSICS_H

#include "particle.h"
#include "magnetic_field.h"
#include "tokamak_geometry.h"
#include <vector>
#include <random>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

class PlasmaPhysics {
private:
    MagneticField& magneticField;
    TokamakGeometry& geometry;

    // Simulation parameters
    float timeScale;
    float plasmaTemperature;
    float particleDensity;
    float velocityScale;

    // Fusion tracking
    float fusionProbability;
    float fusionBoost;
    float maxFusionFractionPerStep;

    // Confinement / stability
    float confinementStrength;
    float coreAttractionStrength;
    float driftOmega;
    float wallLossProbability;
    bool enableCoulomb;
    std::mt19937 rng;

public:
    PlasmaPhysics(MagneticField& field, TokamakGeometry& geom) :
        magneticField(field),
        geometry(geom),
        timeScale(1e-2f),
        plasmaTemperature(1.0e9f),
        particleDensity(1e20f),
        velocityScale(1e-7f),
        fusionProbability(0.0f),
        fusionBoost(1.0e6f),
        maxFusionFractionPerStep(0.02f),
        confinementStrength(50.0f),
        coreAttractionStrength(8.0f),
        driftOmega(2.5f),
        wallLossProbability(0.0f),
        enableCoulomb(false),
        rng(std::random_device{}())
    {}

    // Getters and setters
    float getTimeScale() const { return timeScale; }
    void setTimeScale(float v) { timeScale = v; }
    float getPlasmaTemperature() const { return plasmaTemperature; }
    void setPlasmaTemperature(float v) { plasmaTemperature = v; }
    float getParticleDensity() const { return particleDensity; }
    void setParticleDensity(float v) { particleDensity = v; }
    float getVelocityScale() const { return velocityScale; }
    void setVelocityScale(float v) { velocityScale = v; }
    float getFusionBoost() const { return fusionBoost; }
    void setFusionBoost(float v) { fusionBoost = v; }
    float getMaxFusionFractionPerStep() const { return maxFusionFractionPerStep; }
    void setMaxFusionFractionPerStep(float v) { maxFusionFractionPerStep = v; }
    float getConfinementStrength() const { return confinementStrength; }
    void setConfinementStrength(float v) { confinementStrength = v; }
    float getCoreAttractionStrength() const { return coreAttractionStrength; }
    void setCoreAttractionStrength(float v) { coreAttractionStrength = v; }
    float getDriftOmega() const { return driftOmega; }
    void setDriftOmega(float v) { driftOmega = v; }
    float getWallLossProbability() const { return wallLossProbability; }
    void setWallLossProbability(float v) { wallLossProbability = v; }
    bool getEnableCoulomb() const { return enableCoulomb; }
    void setEnableCoulomb(bool v) { enableCoulomb = v; }

    void updateParticles(std::vector<Particle>& particles, float dt);
    void applyMagneticForce3D(Particle& p, float scaledDt, float realDt);
    void applyCoulombForce(Particle& p1, Particle& p2, float dt);
    bool attemptFusion(Particle& p1, Particle& p2,
                      std::vector<Particle>& newParticles,
                      float dt, bool force);
    void checkBoundaryCollision3D(Particle& p, float dt);
    float getThermalVelocity(float mass) const;
    std::vector<Particle> createThermalPlasma(int numDeuterium, int numTritium);

    /**
     * Inject new fuel particles (D + T) to sustain the reaction.
     * Places them randomly on the torus centerline with radial scatter.
     */
    void injectFuel(std::vector<Particle>& particles, int numD, int numT);
};

// ==================== IMPLEMENTATION ====================

inline void PlasmaPhysics::updateParticles(std::vector<Particle>& particles, float dt)
{
    float scaledDt = dt * timeScale;
    std::vector<Particle> newParticles;

    // Indices for volumetric fusion
    std::vector<size_t> deuteriumIdx;
    std::vector<size_t> tritiumIdx;
    deuteriumIdx.reserve(particles.size());
    tritiumIdx.reserve(particles.size());

    for (size_t i = 0; i < particles.size(); ++i) {
        if (!particles[i].active) continue;

        if (particles[i].type == Particle::DEUTERIUM) deuteriumIdx.push_back(i);
        else if (particles[i].type == Particle::TRITIUM) tritiumIdx.push_back(i);

        // Magnetic confinement (3D)
        applyMagneticForce3D(particles[i], scaledDt, dt);

        // Optional Coulomb interactions
        if (enableCoulomb) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                if (!particles[j].active) continue;
                applyCoulombForce(particles[i], particles[j], scaledDt);
            }
        }

        // Update position (3D)
        particles[i].x += particles[i].vx * scaledDt;
        particles[i].y += particles[i].vy * scaledDt;
        particles[i].z += particles[i].vz * scaledDt;

        // NaN safety
        if (!std::isfinite(particles[i].x) || !std::isfinite(particles[i].y) ||
            !std::isfinite(particles[i].z) ||
            !std::isfinite(particles[i].vx) || !std::isfinite(particles[i].vy) ||
            !std::isfinite(particles[i].vz)) {
            // Reset to a random point on the torus
            float phi = 2.0f * M_PI * (rng() % 10000) / 10000.0f;
            particles[i].x = geometry.torusMajorR * std::cos(phi);
            particles[i].y = 0.0f;
            particles[i].z = geometry.torusMajorR * std::sin(phi);
            particles[i].vx = 0.0f;
            particles[i].vy = 0.0f;
            particles[i].vz = 0.0f;
        }

        // Kinetic energy
        particles[i].kineticEnergy = 0.5f * particles[i].mass *
            (particles[i].vx * particles[i].vx +
             particles[i].vy * particles[i].vy +
             particles[i].vz * particles[i].vz);

        // Boundary collision (3D torus confinement)
        checkBoundaryCollision3D(particles[i], scaledDt);
    }

    // ================== VOLUMETRIC D-T FUSION (TOY MODEL) ==================
    const int ND = (int)deuteriumIdx.size();
    const int NT = (int)tritiumIdx.size();
    const int maxPairs = (ND < NT) ? ND : NT;
    if (maxPairs > 0) {
        // Approximate plasma volume (torus): V = 2π²Rr²
        float R = geometry.torusMajorR;
        float r = geometry.torusMinorR;
        float volume = 2.0f * M_PI * M_PI * R * r * r;
        if (volume < 1e-8f) volume = 1e-8f;

        float nD = (float)ND / volume;
        float nT = (float)NT / volume;

        float T_keV = plasmaTemperature * PhysicsConstants::BOLTZMANN_CONSTANT /
                     (1.0e3f * PhysicsConstants::ELEMENTARY_CHARGE);
        if (T_keV < 1e-6f) T_keV = 1e-6f;
        float reactivity = 1e-6f * std::sqrt(T_keV);

        float expectedFusions = reactivity * nD * nT * volume * dt;
        expectedFusions *= fusionBoost;

        if (expectedFusions > (float)maxPairs) expectedFusions = (float)maxPairs;
        if (expectedFusions < 0.0f) expectedFusions = 0.0f;

        int numFusions = (int)expectedFusions;
        float remainder = expectedFusions - (float)numFusions;
        std::uniform_real_distribution<float> u01(0.0f, 1.0f);
        if (u01(rng) < remainder) numFusions++;

        if (numFusions > maxPairs) numFusions = maxPairs;
        int maxThisStep = (int)((float)maxPairs * maxFusionFractionPerStep);
        if (maxThisStep < 0) maxThisStep = 0;
        if (maxThisStep > maxPairs) maxThisStep = maxPairs;
        if (numFusions > maxThisStep) numFusions = maxThisStep;

        if (numFusions > 0) {
            std::uniform_int_distribution<int> d_pick(0, ND - 1);
            std::uniform_int_distribution<int> t_pick(0, NT - 1);

            for (int k = 0; k < numFusions; ++k) {
                size_t id = deuteriumIdx[(size_t)d_pick(rng)];
                size_t it = tritiumIdx[(size_t)t_pick(rng)];
                if (!particles[id].active || !particles[it].active) continue;
                attemptFusion(particles[id], particles[it], newParticles, scaledDt, true);
            }
        }
    }

    particles.insert(particles.end(), newParticles.begin(), newParticles.end());
}

inline void PlasmaPhysics::applyMagneticForce3D(Particle& p, float scaledDt, float realDt)
{
    if (std::abs(p.charge) < 1e-30f) return;

    // Get 3D magnetic field at particle position
    float Bx, By, Bz;
    magneticField.getTotalField(p.x, p.y, p.z, Bx, By, Bz);

    // Lorentz force: F = q(v × B)
    float Fx, Fy, Fz;
    calculateLorentzForce(p.vx, p.vy, p.vz, Bx, By, Bz, p.charge, Fx, Fy, Fz);

    // Mirror force
    float Fmx, Fmy, Fmz;
    calculateMirrorForce3D(p.x, p.y, p.z, p.vx, p.vy, p.vz, magneticField, p.mass, Fmx, Fmy, Fmz);
    Fx += Fmx;
    Fy += Fmy;
    Fz += Fmz;

    const float forceScale = 1e-6f;
    Fx *= forceScale;
    Fy *= forceScale;
    Fz *= forceScale;

    // Acceleration
    p.vx += (Fx / p.mass) * scaledDt;
    p.vy += (Fy / p.mass) * scaledDt;
    p.vz += (Fz / p.mass) * scaledDt;

    // Core attraction: pull toward the nearest centerline point
    float cx, cy, cz;
    geometry.projectToCenterline(p.x, p.y, p.z, cx, cy, cz);
    float dx = p.x - cx;
    float dy = p.y - cy;
    float dz = p.z - cz;
    float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (dist > 1e-8f) {
        float pull = coreAttractionStrength / (dist + 0.01f);
        p.vx += (-pull * dx) * scaledDt;
        p.vy += (-pull * dy) * scaledDt;
        p.vz += (-pull * dz) * scaledDt;
    }

    // Toroidal drift: rotate particles around the torus ring (in XZ plane)
    float R = std::sqrt(p.x * p.x + p.z * p.z);
    if (R > 1e-6f) {
        // Tangential direction in XZ: (-z/R, 0, x/R)
        float tx = -p.z / R;
        float tz = p.x / R;
        p.vx += driftOmega * tx * scaledDt;
        p.vz += driftOmega * tz * scaledDt;
    }
}

inline void PlasmaPhysics::applyCoulombForce(Particle& p1, Particle& p2, float dt)
{
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float dz = p2.z - p1.z;
    float r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1e-6f) r = 1e-6f;

    float forceMagnitude = PhysicsConstants::COULOMB_CONSTANT *
                           p1.charge * p2.charge / (r * r);

    float fx = forceMagnitude * dx / r;
    float fy = forceMagnitude * dy / r;
    float fz = forceMagnitude * dz / r;

    float debyeLength = 7.43e2f * std::sqrt(plasmaTemperature / particleDensity);
    float screeningFactor = std::exp(-r / debyeLength);
    fx *= screeningFactor;
    fy *= screeningFactor;
    fz *= screeningFactor;

    const float forceScale = 1e-6f;
    fx *= forceScale;
    fy *= forceScale;
    fz *= forceScale;

    p1.vx += (fx / p1.mass) * dt;
    p1.vy += (fy / p1.mass) * dt;
    p1.vz += (fz / p1.mass) * dt;
    p2.vx += (-fx / p2.mass) * dt;
    p2.vy += (-fy / p2.mass) * dt;
    p2.vz += (-fz / p2.mass) * dt;
}

inline bool PlasmaPhysics::attemptFusion(Particle& p1, Particle& p2,
    std::vector<Particle>& newParticles, float dt, bool force)
{
    float vrel_x = p1.vx - p2.vx;
    float vrel_y = p1.vy - p2.vy;
    float vrel_z = p1.vz - p2.vz;
    float vrel = std::sqrt(vrel_x * vrel_x + vrel_y * vrel_y + vrel_z * vrel_z);

    float reducedMass = (p1.mass * p2.mass) / (p1.mass + p2.mass);
    float E_cm = 0.5f * reducedMass * vrel * vrel;

    if (!force) {
        if (E_cm < PhysicsConstants::FUSION_THRESHOLD_ENERGY) return false;
    }

    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float dz = p2.z - p1.z;
    float distance = std::sqrt(dx * dx + dy * dy + dz * dz);

    float crossSection = PhysicsConstants::FUSION_CROSS_SECTION *
                        (E_cm / PhysicsConstants::FUSION_THRESHOLD_ENERGY);

    if (!force) {
        const float fusionInteractionDistance = 0.03f;
        if (distance > fusionInteractionDistance) return false;
        float fusionChance = crossSection * particleDensity * vrel * dt;
        fusionChance *= fusionBoost;
        if (fusionChance < 0.0f) fusionChance = 0.0f;
        if (fusionChance > 1.0f) fusionChance = 1.0f;
        std::uniform_real_distribution<float> dist01(0.0f, 1.0f);
        if (dist01(rng) > fusionChance) return false;
    }

    // === FUSION: D + T → He⁴ + n + 17.6 MeV ===
    float cm_x = (p1.mass * p1.x + p2.mass * p2.x) / (p1.mass + p2.mass);
    float cm_y = (p1.mass * p1.y + p2.mass * p2.y) / (p1.mass + p2.mass);
    float cm_z = (p1.mass * p1.z + p2.mass * p2.z) / (p1.mass + p2.mass);
    float cm_vx = (p1.mass * p1.vx + p2.mass * p2.vx) / (p1.mass + p2.mass);
    float cm_vy = (p1.mass * p1.vy + p2.mass * p2.vy) / (p1.mass + p2.mass);
    float cm_vz = (p1.mass * p1.vz + p2.mass * p2.vz) / (p1.mass + p2.mass);

    float E_alpha = 3.5e6f * PhysicsConstants::ELEMENTARY_CHARGE;
    float E_neutron = 14.1e6f * PhysicsConstants::ELEMENTARY_CHARGE;

    // Random 3D direction for product emission
    std::uniform_real_distribution<float> angle_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> cos_dist(-1.0f, 1.0f);
    float phi = angle_dist(rng);
    float cosTheta = cos_dist(rng);
    float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);

    float v_alpha = std::sqrt(2.0f * E_alpha / PhysicsConstants::HELIUM_MASS);
    float v_neutron = std::sqrt(2.0f * E_neutron / PhysicsConstants::NEUTRON_MASS);

    float dirx = sinTheta * std::cos(phi);
    float diry = sinTheta * std::sin(phi);
    float dirz = cosTheta;

    float vx_he = (cm_vx + v_alpha * dirx) * velocityScale;
    float vy_he = (cm_vy + v_alpha * diry) * velocityScale;
    float vz_he = (cm_vz + v_alpha * dirz) * velocityScale;
    float vx_n = (cm_vx - v_neutron * dirx) * velocityScale;
    float vy_n = (cm_vy - v_neutron * diry) * velocityScale;
    float vz_n = (cm_vz - v_neutron * dirz) * velocityScale;

    Particle helium = createParticle(Particle::HELIUM, cm_x, cm_y, vx_he, vy_he, cm_z, vz_he);
    Particle neutron = createParticle(Particle::NEUTRON, cm_x, cm_y, vx_n, vy_n, cm_z, vz_n);

    newParticles.push_back(helium);
    newParticles.push_back(neutron);

    p1.active = false;
    p2.active = false;

    return true;
}

inline void PlasmaPhysics::checkBoundaryCollision3D(Particle& p, float dt)
{
    float sdf = geometry.torusSDF(p.x, p.y, p.z);

    if (sdf > 0.0f) {
        // Particle outside torus — push it back
        float nx, ny, nz;
        geometry.torusNormal(p.x, p.y, p.z, nx, ny, nz);

        // Inward acceleration
        float pushStrength = confinementStrength * sdf;
        p.vx -= pushStrength * nx * dt;
        p.vy -= pushStrength * ny * dt;
        p.vz -= pushStrength * nz * dt;

        // Hard constraint: project back onto torus boundary
        float edgeBuffer = 0.01f;
        p.x -= (sdf + edgeBuffer) * nx * 1.05f;
        p.y -= (sdf + edgeBuffer) * ny * 1.05f;
        p.z -= (sdf + edgeBuffer) * nz * 1.05f;

        // Remove outward velocity component
        float vdotn = p.vx * nx + p.vy * ny + p.vz * nz;
        if (vdotn > 0.0f) {
            p.vx -= vdotn * nx;
            p.vy -= vdotn * ny;
            p.vz -= vdotn * nz;
        }

        // Wall losses
        if (wallLossProbability > 0.0f) {
            std::uniform_real_distribution<float> u01(0.0f, 1.0f);
            if (u01(rng) < wallLossProbability) {
                p.active = false;
            }
        }
    } else if (sdf > -0.02f) {
        // Close to boundary — gentle inward push
        float nx, ny, nz;
        geometry.torusNormal(p.x, p.y, p.z, nx, ny, nz);
        float penetration = sdf + 0.02f;
        p.vx -= confinementStrength * penetration * nx * dt;
        p.vy -= confinementStrength * penetration * ny * dt;
        p.vz -= confinementStrength * penetration * nz * dt;

        float vdotn = p.vx * nx + p.vy * ny + p.vz * nz;
        if (vdotn > 0.0f) {
            p.vx -= vdotn * nx;
            p.vy -= vdotn * ny;
            p.vz -= vdotn * nz;
        }
    }
}

inline float PlasmaPhysics::getThermalVelocity(float mass) const
{
    return std::sqrt(3.0f * PhysicsConstants::BOLTZMANN_CONSTANT *
                     plasmaTemperature / mass);
}

inline std::vector<Particle> PlasmaPhysics::createThermalPlasma(
    int numDeuterium, int numTritium)
{
    std::vector<Particle> particles;
    std::uniform_real_distribution<float> phi_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> theta_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> r_dist(0.0f, 1.0f);
    std::normal_distribution<float> vel_d(0.0f, getThermalVelocity(PhysicsConstants::DEUTERIUM_MASS) * velocityScale);
    std::normal_distribution<float> vel_t(0.0f, getThermalVelocity(PhysicsConstants::TRITIUM_MASS) * velocityScale);

    float R = geometry.torusMajorR;
    float rr = geometry.torusMinorR;

    // Create deuterium particles distributed inside the torus
    for (int i = 0; i < numDeuterium; ++i) {
        float phi = phi_dist(rng);   // toroidal angle
        float theta = theta_dist(rng); // poloidal angle
        float rFrac = std::sqrt(r_dist(rng)) * rr * 0.85f; // concentrate toward center

        float x = (R + rFrac * std::cos(theta)) * std::cos(phi);
        float y = rFrac * std::sin(theta);
        float z = (R + rFrac * std::cos(theta)) * std::sin(phi);

        float vx = vel_d(rng);
        float vy = vel_d(rng);
        float vz = vel_d(rng);

        particles.push_back(createParticle(Particle::DEUTERIUM, x, y, vx, vy, z, vz));
    }

    // Create tritium particles
    for (int i = 0; i < numTritium; ++i) {
        float phi = phi_dist(rng);
        float theta = theta_dist(rng);
        float rFrac = std::sqrt(r_dist(rng)) * rr * 0.85f;

        float x = (R + rFrac * std::cos(theta)) * std::cos(phi);
        float y = rFrac * std::sin(theta);
        float z = (R + rFrac * std::cos(theta)) * std::sin(phi);

        float vx = vel_t(rng);
        float vy = vel_t(rng);
        float vz = vel_t(rng);

        particles.push_back(createParticle(Particle::TRITIUM, x, y, vx, vy, z, vz));
    }

    return particles;
}

inline void PlasmaPhysics::injectFuel(std::vector<Particle>& particles, int numD, int numT)
{
    std::uniform_real_distribution<float> phi_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> theta_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> r_dist(0.0f, 1.0f);
    std::normal_distribution<float> vel_d(0.0f, getThermalVelocity(PhysicsConstants::DEUTERIUM_MASS) * velocityScale);
    std::normal_distribution<float> vel_t(0.0f, getThermalVelocity(PhysicsConstants::TRITIUM_MASS) * velocityScale);

    float R = geometry.torusMajorR;
    float rr = geometry.torusMinorR;

    for (int i = 0; i < numD; ++i) {
        float phi = phi_dist(rng);
        float theta = theta_dist(rng);
        float rFrac = std::sqrt(r_dist(rng)) * rr * 0.7f;

        float x = (R + rFrac * std::cos(theta)) * std::cos(phi);
        float y = rFrac * std::sin(theta);
        float z = (R + rFrac * std::cos(theta)) * std::sin(phi);

        particles.push_back(createParticle(Particle::DEUTERIUM, x, y, vel_d(rng), vel_d(rng), z, vel_d(rng)));
    }

    for (int i = 0; i < numT; ++i) {
        float phi = phi_dist(rng);
        float theta = theta_dist(rng);
        float rFrac = std::sqrt(r_dist(rng)) * rr * 0.7f;

        float x = (R + rFrac * std::cos(theta)) * std::cos(phi);
        float y = rFrac * std::sin(theta);
        float z = (R + rFrac * std::cos(theta)) * std::sin(phi);

        particles.push_back(createParticle(Particle::TRITIUM, x, y, vel_t(rng), vel_t(rng), z, vel_t(rng)));
    }
}

#endif // PLASMA_PHYSICS_H
