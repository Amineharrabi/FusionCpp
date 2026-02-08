#ifndef PLASMA_PHYSICS_H
#define PLASMA_PHYSICS_H

#include "particle.h"
#include "magnetic_field.h"
#include "tokamak_geometry.h"
#include <vector>
#include <random>
#include <cmath>
#define M_PI 3.141592


class PlasmaPhysics {
private:
MagneticField& magneticField;
TokamakGeometry& geometry;

    // Simulation parameters
    float timeScale;        // Time scaling factor for stability
    float plasmaTemperature; // Temperature in Kelvin (typical: 10⁸ K = 10 keV)
    float particleDensity;   // Particles per m³
    float velocityScale;

    // Collision and fusion tracking
    float fusionProbability;
    float fusionBoost;
    float maxFusionFractionPerStep;

    // Confinement / stability knobs (toy model)
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
plasmaTemperature(1.0e9f), // Hotter plasma so D-T fusion is reachable in this toy model
particleDensity(1e20f), // 10²⁰ particles/m³
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

    /**
     * Update all particles for one time step
     * Applies forces: magnetic confinement, Coulomb repulsion, fusion
     */
    void updateParticles(std::vector<Particle>& particles, float dt);

    /**
     * Apply magnetic confinement force to a single particle
     * Uses Lorentz force and mirror force
     */
    void applyMagneticForce(Particle& p, float scaledDt, float realDt);

    /**
     * Calculate Coulomb force between two charged particles
     * F = k·q₁·q₂·r̂ / r²
     */
    void applyCoulombForce(Particle& p1, Particle& p2, float dt);

    /**
     * Check and execute fusion reaction between particles
     * D + T → He⁴ + n + 17.6 MeV
     */
    bool attemptFusion(Particle& p1, Particle& p2,
                      std::vector<Particle>& newParticles,
                      float dt,
                      bool force);

    /**
     * Handle particle collision with plasma boundary
     * Particles hitting the wall lose energy or reflect
     */
    void checkBoundaryCollision(Particle& p, float dt);

    /**
     * Calculate thermal velocity from temperature
     * v_thermal = sqrt(3kT/m)
     */
    float getThermalVelocity(float mass) const;

    /**
     * Initialize particles with thermal distribution
     */
    std::vector<Particle> createThermalPlasma(int numDeuterium, int numTritium);

};

inline void PlasmaPhysics::updateParticles(std::vector<Particle>& particles, float dt)
{
float scaledDt = dt * timeScale;
std::vector<Particle> newParticles;

    // Count active fuel particles for volumetric fusion model
    std::vector<size_t> deuteriumIdx;
    std::vector<size_t> tritiumIdx;
    deuteriumIdx.reserve(particles.size());
    tritiumIdx.reserve(particles.size());

    // Apply forces to each active particle
    for (size_t i = 0; i < particles.size(); ++i) {
        if (!particles[i].active) continue;

        if (particles[i].type == Particle::DEUTERIUM) deuteriumIdx.push_back(i);
        else if (particles[i].type == Particle::TRITIUM) tritiumIdx.push_back(i);

        // Magnetic confinement
        applyMagneticForce(particles[i], scaledDt, dt);

        // Optional Coulomb interactions (very expensive O(N^2); disabled by default)
        if (enableCoulomb) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                if (!particles[j].active) continue;
                applyCoulombForce(particles[i], particles[j], scaledDt);
            }
        }

        // Update position
        particles[i].x += particles[i].vx * scaledDt;
        particles[i].y += particles[i].vy * scaledDt;

        // Numerical safety: keep simulation from producing NaNs/Infs that make particles disappear
        if (!std::isfinite(particles[i].x) || !std::isfinite(particles[i].y) ||
            !std::isfinite(particles[i].vx) || !std::isfinite(particles[i].vy)) {
            particles[i].x = 0.0f;
            particles[i].y = 0.0f;
            particles[i].vx = 0.0f;
            particles[i].vy = 0.0f;
        }

        // Update kinetic energy
        particles[i].kineticEnergy = 0.5f * particles[i].mass *
            (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);

        // Check boundary collisions
        checkBoundaryCollision(particles[i], scaledDt);
    }

    // ================== VOLUMETRIC D-T FUSION (TOY MODEL) ==================
    // Use a temperature-driven reactivity and fuse random D/T pairs, rather than relying
    // on rare close encounters in 2D.
    const int ND = (int)deuteriumIdx.size();
    const int NT = (int)tritiumIdx.size();
    const int maxPairs = (ND < NT) ? ND : NT;
    if (maxPairs > 0) {
        // Estimate plasma cross-section area from the drawn boundary polygon
        float area = 0.0f;
        const size_t n = geometry.plasmaVertices.size() / 2;
        if (n >= 3) {
            for (size_t i = 0; i < n; ++i) {
                size_t j = (i + 1) % n;
                float xi = geometry.plasmaVertices[i * 2 + 0];
                float yi = geometry.plasmaVertices[i * 2 + 1];
                float xj = geometry.plasmaVertices[j * 2 + 0];
                float yj = geometry.plasmaVertices[j * 2 + 1];
                area += (xi * yj - xj * yi);
            }
            area = std::abs(area) * 0.5f;
        }
        if (area > 1e-8f) {
            // 2D "number densities" (particles per unit area). This is a visualization model.
            float nD_2d = (float)ND / area;
            float nT_2d = (float)NT / area;

            // Very rough DT reactivity proxy increasing with temperature.
            // Keep it simple for tuning: <σv> ~ C * sqrt(T_keV)
            float T_keV = plasmaTemperature * PhysicsConstants::BOLTZMANN_CONSTANT / (1.0e3f * PhysicsConstants::ELEMENTARY_CHARGE);
            if (T_keV < 1e-6f) T_keV = 1e-6f;
            float reactivity = 1e-6f * std::sqrt(T_keV); // visualization units

            // IMPORTANT: fusion rate is integrated in real (frame) dt, not visualization-scaled dt,
            // otherwise changing timeScale would incorrectly change fusion frequency.
            float expectedFusions = reactivity * nD_2d * nT_2d * area * dt;
            expectedFusions *= fusionBoost;

            // Clamp expected to avoid pathological spikes
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

                    // Fuse the two selected fuel ions
                    attemptFusion(particles[id], particles[it], newParticles, scaledDt, true);
                }
            }
        }
    }

    // Add newly created particles (from fusion)
    particles.insert(particles.end(), newParticles.begin(), newParticles.end());

}

inline void PlasmaPhysics::applyMagneticForce(Particle& p, float scaledDt, float realDt)
{
if (std::abs(p.charge) < 1e-30f) return; // Neutral particles not affected

    // Get magnetic field at particle position
    float Bx, By, Bz;
    magneticField.getTotalField(p.x, p.y, Bx, By, Bz);

    // Calculate Lorentz force: F = q(v × B)
    float Fx, Fy, Fz;
    calculateLorentzForce(p.vx, p.vy, 0.0f, Bx, By, Bz, p.charge, Fx, Fy, Fz);

    // Add magnetic mirror force (for confinement in strong field regions)
    float Fx_mirror, Fy_mirror;
    calculateMirrorForce(p.x, p.y, p.vx, p.vy, magneticField, p.mass, Fx_mirror, Fy_mirror);

    Fx += Fx_mirror;
    Fy += Fy_mirror;

    // This simulation uses screen-space geometry (order ~1) but SI charges/masses, which can
    // produce enormous accelerations. Scale forces down for a stable visualization.
    const float forceScale = 1e-6f;
    Fx *= forceScale;
    Fy *= forceScale;

    // Apply acceleration: a = F/m
    float ax = Fx / p.mass;
    float ay = Fy / p.mass;

    // Update velocity
    p.vx += ax * scaledDt;
    p.vy += ay * scaledDt;

    // Core-localizing terms for 2D visualization
    float rx = p.x;
    float ry = p.y;
    float r2 = rx * rx + ry * ry;
    if (r2 > 1e-12f) {
        // Gentle inward pull toward core
        p.vx += (-coreAttractionStrength * rx) * scaledDt;
        p.vy += (-coreAttractionStrength * ry) * scaledDt;
    }

    // Rotational drift around the core (circulation)
    p.vx += (-driftOmega * p.y) * scaledDt;
    p.vy += ( driftOmega * p.x) * scaledDt;

}

inline void PlasmaPhysics::applyCoulombForce(Particle& p1, Particle& p2, float dt)
{
// Calculate distance between particles
float dx = p2.x - p1.x;
float dy = p2.y - p1.y;
float r = std::sqrt(dx*dx + dy*dy);

    // Avoid singularity at r=0
    if (r < 1e-6f) r = 1e-6f;

    // Coulomb force: F = k·q₁·q₂ / r²
    // k = 8.988×10⁹ N⋅m²/C²
    float forceMagnitude = PhysicsConstants::COULOMB_CONSTANT *
                           p1.charge * p2.charge / (r * r);

    // Force direction (unit vector)
    float fx = forceMagnitude * dx / r;
    float fy = forceMagnitude * dy / r;

    // Apply screening for plasma (Debye shielding)
    // Force reduced at distances > Debye length
    float debyeLength = 7.43e2f * std::sqrt(plasmaTemperature / particleDensity);
    float screeningFactor = std::exp(-r / debyeLength);
    fx *= screeningFactor;
    fy *= screeningFactor;

    // Scale down Coulomb interaction for stability/visualization
    const float forceScale = 1e-6f;
    fx *= forceScale;
    fy *= forceScale;

    // Apply forces (Newton's third law)
    float ax1 = fx / p1.mass;
    float ay1 = fy / p1.mass;
    float ax2 = -fx / p2.mass;
    float ay2 = -fy / p2.mass;

    p1.vx += ax1 * dt;
    p1.vy += ay1 * dt;
    p2.vx += ax2 * dt;
    p2.vy += ay2 * dt;

}

inline bool PlasmaPhysics::attemptFusion(Particle& p1, Particle& p2,
std::vector<Particle>& newParticles,
float dt,
bool force)
{
// Calculate relative velocity
float vrel_x = p1.vx - p2.vx;
float vrel_y = p1.vy - p2.vy;
float vrel = std::sqrt(vrel_x*vrel_x + vrel_y*vrel_y);

    // Calculate center of mass energy
    float reducedMass = (p1.mass * p2.mass) / (p1.mass + p2.mass);
    float E_cm = 0.5f * reducedMass * vrel * vrel;

    // Check if energy exceeds fusion threshold (~60 keV for D-T)
    if (!force) {
        if (E_cm < PhysicsConstants::FUSION_THRESHOLD_ENERGY) {
            return false;
        }
    }

    // Calculate distance
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float distance = std::sqrt(dx*dx + dy*dy);

    // Fusion probability based on cross-section σ(E)
    // P ≈ σ(E) · n · v · dt
    float crossSection = PhysicsConstants::FUSION_CROSS_SECTION *
                        (E_cm / PhysicsConstants::FUSION_THRESHOLD_ENERGY);

    if (!force) {
        // Visualization-friendly interaction distance (this is not physical; it's a toy model)
        const float fusionInteractionDistance = 0.01f;
        if (distance > fusionInteractionDistance) return false;

        float fusionChance = crossSection * particleDensity * vrel * dt;
        fusionChance *= fusionBoost;
        if (fusionChance < 0.0f) fusionChance = 0.0f;
        if (fusionChance > 1.0f) fusionChance = 1.0f;

        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        float randomValue = dist(rng);
        if (randomValue > fusionChance) return false;
    }

    // === FUSION OCCURRED: D + T → He⁴ + n + 17.6 MeV ===

    // Center of mass
    float cm_x = (p1.mass * p1.x + p2.mass * p2.x) / (p1.mass + p2.mass);
    float cm_y = (p1.mass * p1.y + p2.mass * p2.y) / (p1.mass + p2.mass);
    float cm_vx = (p1.mass * p1.vx + p2.mass * p2.vx) / (p1.mass + p2.mass);
    float cm_vy = (p1.mass * p1.vy + p2.mass * p2.vy) / (p1.mass + p2.mass);

    // Energy release: 17.6 MeV = 2.82×10⁻¹² J
    float Q_fusion = 2.82e-12f;

    // Products: He-4 (alpha particle) gets ~3.5 MeV, neutron gets ~14.1 MeV
    float E_alpha = 3.5e6f * PhysicsConstants::ELEMENTARY_CHARGE; // Convert eV to J
    float E_neutron = 14.1e6f * PhysicsConstants::ELEMENTARY_CHARGE;

    // Calculate product velocities (isotropic emission in CM frame)
    std::uniform_real_distribution<float> angle_dist(0.0f, 2.0f * M_PI);
    float theta = angle_dist(rng);

    // Alpha particle
    float v_alpha = std::sqrt(2.0f * E_alpha / PhysicsConstants::HELIUM_MASS);
    float v_neutron = std::sqrt(2.0f * E_neutron / PhysicsConstants::NEUTRON_MASS);

    float vx_he = (cm_vx + v_alpha * std::cos(theta)) * velocityScale;
    float vy_he = (cm_vy + v_alpha * std::sin(theta)) * velocityScale;
    float vx_n = (cm_vx - v_neutron * std::cos(theta)) * velocityScale;
    float vy_n = (cm_vy - v_neutron * std::sin(theta)) * velocityScale;

    Particle helium = createParticle(Particle::HELIUM, cm_x, cm_y, vx_he, vy_he);
    helium.radius = 0.005f; // Slightly larger for visibility

    // Neutron (opposite direction for momentum conservation)
    Particle neutron = createParticle(Particle::NEUTRON, cm_x, cm_y, vx_n, vy_n);

    newParticles.push_back(helium);
    newParticles.push_back(neutron);

    // Deactivate reactant particles
    p1.active = false;
    p2.active = false;

    return true;

}

inline void PlasmaPhysics::checkBoundaryCollision(Particle& p, float dt)
{
float distFromEdge = geometry.distanceFromPlasmaEdge(p.x, p.y);

    // Keep particles away from the separatrix (visual buffer inside the plasma)
    const float edgeBuffer = 0.015f;

    // Compute boundary normal from signed distance gradient (points outward)
    float dx = 0.001f;
    float d1 = geometry.distanceFromPlasmaEdge(p.x - dx, p.y);
    float d2 = geometry.distanceFromPlasmaEdge(p.x + dx, p.y);
    float nx = (d2 - d1) / (2.0f * dx);

    d1 = geometry.distanceFromPlasmaEdge(p.x, p.y - dx);
    d2 = geometry.distanceFromPlasmaEdge(p.x, p.y + dx);
    float ny = (d2 - d1) / (2.0f * dx);

    float norm = std::sqrt(nx * nx + ny * ny) + 1e-10f;
    nx /= norm;
    ny /= norm;

    // Particle outside plasma boundary
    if (distFromEdge > 0.0f) {
        // Inward acceleration proportional to penetration distance
        float ax = -confinementStrength * distFromEdge * nx;
        float ay = -confinementStrength * distFromEdge * ny;
        p.vx += ax * dt;
        p.vy += ay * dt;

        // Hard constraint for confinement: project particle back onto the boundary and
        // remove outward normal velocity component (no bounce, no wall heating here).
        p.x -= (distFromEdge + edgeBuffer) * nx * 1.05f;
        p.y -= (distFromEdge + edgeBuffer) * ny * 1.05f;

        float vdotn = p.vx * nx + p.vy * ny;
        if (vdotn > 0.0f) {
            p.vx -= vdotn * nx;
            p.vy -= vdotn * ny;
        }

        // Optional: wall losses (particles hitting vessel are lost)
        if (wallLossProbability > 0.0f) {
            std::uniform_real_distribution<float> u01(0.0f, 1.0f);
            float lossChance = wallLossProbability;
            if (lossChance > 1.0f) lossChance = 1.0f;
            if (u01(rng) < lossChance) {
                p.active = false;
            }
        }
    } else if (distFromEdge > -edgeBuffer) {
        // Particle is inside but too close to the boundary: apply a gentle inward push
        // so particles "float" in the plasma volume instead of sticking to the edge.
        float penetration = distFromEdge + edgeBuffer; // in [0, edgeBuffer]
        float ax = -confinementStrength * penetration * nx;
        float ay = -confinementStrength * penetration * ny;
        p.vx += ax * dt;
        p.vy += ay * dt;

        // Remove outward normal component so the edge acts like a magnetic surface
        float vdotn = p.vx * nx + p.vy * ny;
        if (vdotn > 0.0f) {
            p.vx -= vdotn * nx;
            p.vy -= vdotn * ny;
        }
    }

}

inline float PlasmaPhysics::getThermalVelocity(float mass) const
{
// v_thermal = sqrt(3kT/m)
return std::sqrt(3.0f * PhysicsConstants::BOLTZMANN_CONSTANT *
plasmaTemperature / mass);
}

inline std::vector<Particle> PlasmaPhysics::createThermalPlasma(
int numDeuterium, int numTritium)
{
std::vector<Particle> particles;
std::uniform_real_distribution<float> pos_dist(-0.2f, 0.2f);
std::normal_distribution<float> vel_dist_d(0.0f, getThermalVelocity(PhysicsConstants::DEUTERIUM_MASS) * velocityScale);
std::normal_distribution<float> vel_dist_t(0.0f, getThermalVelocity(PhysicsConstants::TRITIUM_MASS) * velocityScale);

    // Create deuterium particles
    for (int i = 0; i < numDeuterium; ++i) {
        float x, y;
        do {
            x = pos_dist(rng);
            y = pos_dist(rng);
        } while (!geometry.isInsidePlasma(x, y));

        float vx = vel_dist_d(rng);
        float vy = vel_dist_d(rng);

        particles.push_back(createParticle(Particle::DEUTERIUM, x, y, vx, vy));
    }

    // Create tritium particles
    for (int i = 0; i < numTritium; ++i) {
        float x, y;
        do {
            x = pos_dist(rng);
            y = pos_dist(rng);
        } while (!geometry.isInsidePlasma(x, y));

        float vx = vel_dist_t(rng);
        float vy = vel_dist_t(rng);

        particles.push_back(createParticle(Particle::TRITIUM, x, y, vx, vy));
    }

    return particles;

}

#endif // PLASMA_PHYSICS_H
