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

    // Collision and fusion tracking
    float fusionProbability;
    std::mt19937 rng;

public:
PlasmaPhysics(MagneticField& field, TokamakGeometry& geom) :
magneticField(field),
geometry(geom),
timeScale(1e-9f), // Nanosecond scale for stability
plasmaTemperature(1.5e8f), // 150 million Kelvin (~13 keV)
particleDensity(1e20f), // 10²⁰ particles/m³
fusionProbability(0.0f),
rng(std::random_device{}())
{}

    /**
     * Update all particles for one time step
     * Applies forces: magnetic confinement, Coulomb repulsion, fusion
     */
    void updateParticles(std::vector<Particle>& particles, float dt);

    /**
     * Apply magnetic confinement force to a single particle
     * Uses Lorentz force and mirror force
     */
    void applyMagneticForce(Particle& p, float dt);

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
                      std::vector<Particle>& newParticles);

    /**
     * Handle particle collision with plasma boundary
     * Particles hitting the wall lose energy or reflect
     */
    void checkBoundaryCollision(Particle& p);

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

    // Apply forces to each active particle
    for (size_t i = 0; i < particles.size(); ++i) {
        if (!particles[i].active) continue;

        // Magnetic confinement
        applyMagneticForce(particles[i], scaledDt);

        // Coulomb interactions with other particles
        for (size_t j = i + 1; j < particles.size(); ++j) {
            if (!particles[j].active) continue;

            applyCoulombForce(particles[i], particles[j], scaledDt);

            // Check for fusion
            if ((particles[i].type == Particle::DEUTERIUM && particles[j].type == Particle::TRITIUM) ||
                (particles[i].type == Particle::TRITIUM && particles[j].type == Particle::DEUTERIUM)) {

                if (attemptFusion(particles[i], particles[j], newParticles)) {
                    // Fusion occurred, particles consumed
                    continue;
                }
            }
        }

        // Update position
        particles[i].x += particles[i].vx * scaledDt;
        particles[i].y += particles[i].vy * scaledDt;

        // Update kinetic energy
        particles[i].kineticEnergy = 0.5f * particles[i].mass *
            (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);

        // Check boundary collisions
        checkBoundaryCollision(particles[i]);
    }

    // Add newly created particles (from fusion)
    particles.insert(particles.end(), newParticles.begin(), newParticles.end());

}

inline void PlasmaPhysics::applyMagneticForce(Particle& p, float dt)
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

    // Apply acceleration: a = F/m
    float ax = Fx / p.mass;
    float ay = Fy / p.mass;

    // Update velocity
    p.vx += ax * dt;
    p.vy += ay * dt;

    // Apply damping to prevent numerical instability
    float damping = 0.9999f;
    p.vx *= damping;
    p.vy *= damping;

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
std::vector<Particle>& newParticles)
{
// Calculate relative velocity
float vrel_x = p1.vx - p2.vx;
float vrel_y = p1.vy - p2.vy;
float vrel = std::sqrt(vrel_x*vrel_x + vrel_y*vrel_y);

    // Calculate center of mass energy
    float reducedMass = (p1.mass * p2.mass) / (p1.mass + p2.mass);
    float E_cm = 0.5f * reducedMass * vrel * vrel;

    // Check if energy exceeds fusion threshold (~60 keV for D-T)
    if (E_cm < PhysicsConstants::FUSION_THRESHOLD_ENERGY) {
        return false;
    }

    // Calculate distance
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float distance = std::sqrt(dx*dx + dy*dy);

    // Fusion probability based on cross-section σ(E)
    // P = σ·n·v·dt where σ increases with energy
    float crossSection = PhysicsConstants::FUSION_CROSS_SECTION *
                        (E_cm / PhysicsConstants::FUSION_THRESHOLD_ENERGY);

    // Simple distance-based probability
    if (distance > 1e-4f) return false;

    // Random fusion probability
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float randomValue = dist(rng);
    float fusionChance = crossSection * 1e28f; // Scaled for simulation

    if (randomValue > fusionChance) {
        return false;
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
    Particle alpha = createParticle(Particle::HELIUM, cm_x, cm_y,
                                   cm_vx + v_alpha * std::cos(theta),
                                   cm_vy + v_alpha * std::sin(theta));
    alpha.radius = 0.005f; // Slightly larger for visibility

    // Neutron (opposite direction for momentum conservation)
    float v_neutron = std::sqrt(2.0f * E_neutron / PhysicsConstants::NEUTRON_MASS);
    Particle neutron = createParticle(Particle::NEUTRON, cm_x, cm_y,
                                     cm_vx - v_neutron * std::cos(theta),
                                     cm_vy - v_neutron * std::sin(theta));

    newParticles.push_back(alpha);
    newParticles.push_back(neutron);

    // Deactivate reactant particles
    p1.active = false;
    p2.active = false;

    return true;

}

inline void PlasmaPhysics::checkBoundaryCollision(Particle& p)
{
float distFromEdge = geometry.distanceFromPlasmaEdge(p.x, p.y);

    // Particle outside plasma boundary
    if (distFromEdge > 0.0f) {
        // Reflect particle back (elastic collision with wall)
        // Calculate normal to boundary
        float dx = 0.001f;
        float d1 = geometry.distanceFromPlasmaEdge(p.x - dx, p.y);
        float d2 = geometry.distanceFromPlasmaEdge(p.x + dx, p.y);
        float nx = (d2 - d1) / (2.0f * dx);

        d1 = geometry.distanceFromPlasmaEdge(p.x, p.y - dx);
        d2 = geometry.distanceFromPlasmaEdge(p.x, p.y + dx);
        float ny = (d2 - d1) / (2.0f * dx);

        // Normalize
        float norm = std::sqrt(nx*nx + ny*ny) + 1e-10f;
        nx /= norm;
        ny /= norm;

        // Reflect velocity: v' = v - 2(v·n)n
        float vdotn = p.vx * nx + p.vy * ny;
        p.vx -= 2.0f * vdotn * nx;
        p.vy -= 2.0f * vdotn * ny;

        // Energy loss on collision (wall absorption)
        float energyLoss = 0.95f;
        p.vx *= energyLoss;
        p.vy *= energyLoss;

        // Move particle back inside
        p.x -= distFromEdge * nx * 1.1f;
        p.y -= distFromEdge * ny * 1.1f;
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
std::normal_distribution<float> vel_dist_d(0.0f, getThermalVelocity(PhysicsConstants::DEUTERIUM_MASS));
std::normal_distribution<float> vel_dist_t(0.0f, getThermalVelocity(PhysicsConstants::TRITIUM_MASS));

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


#endif // PLASMA_PHYSICS_H
