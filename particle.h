#ifndef PARTICLE_H
#define PARTICLE_H

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif
#include <vector>
#include <cmath>

struct Particle
{
    // Position (in simulation space)
    float x, y, z;

    // Velocity
    float vx, vy, vz;

    // Physical properties
    float mass;   // kg
    float charge; // Coulombs
    float radius; // Visual radius for rendering

    // Visual properties
    float r, g, b, a; // RGBA color

    // Particle type identifier
    enum Type
    {
        DEUTERIUM,
        TRITIUM,
        HELIUM,
        NEUTRON,
        ELECTRON
    } type;

    // Energy for fusion calculations (Joules)
    float kineticEnergy;

    // Flag for particles that have fused
    bool active;
};

// GPU-side particle for SSBO upload (std430 aligned)
struct GPUParticle {
    float px, py, pz, radius;   // posRadius
    float r, g, b, a;           // color
};

// Fusion flash event for GPU rendering
struct FusionFlash {
    float px, py, pz, age;      // posAge  (age 0..1, 0=new, 1=faded)
    float r, g, b, intensity;   // color + intensity
};

// Physical constants
namespace PhysicsConstants
{
    constexpr float ELECTRON_MASS = 9.109e-31f;  // kg
    constexpr float PROTON_MASS = 1.673e-27f;    // kg
    constexpr float DEUTERIUM_MASS = 3.344e-27f; // kg
    constexpr float TRITIUM_MASS = 5.008e-27f;   // kg
    constexpr float HELIUM_MASS = 6.646e-27f;    // kg
    constexpr float NEUTRON_MASS = 1.675e-27f;   // kg

    constexpr float ELEMENTARY_CHARGE = 1.602e-19f;   // Coulombs
    constexpr float VACUUM_PERMITTIVITY = 8.854e-12f; // F/m
    constexpr float COULOMB_CONSTANT = 8.988e9f;      // N·m²/C²
    constexpr float BOLTZMANN_CONSTANT = 1.381e-23f;  // J/K

    // Fusion cross-section parameters (simplified)
    constexpr float FUSION_THRESHOLD_ENERGY = 1.0e-14f; // ~60 keV in Joules
    constexpr float FUSION_CROSS_SECTION = 1.0e-28f;    // m²
}

// Create particle with default properties based on type (3D version)
inline Particle createParticle(Particle::Type type, float x, float y, float vx, float vy,
                               float z = 0.0f, float vz = 0.0f)
{
    Particle p;
    p.x = x;
    p.y = y;
    p.z = z;
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
    p.type = type;
    p.active = true;
    p.radius = 0.02f; // 3D visual radius (larger for visibility in 3D)

    switch (type)
    {
    case Particle::DEUTERIUM:
        p.mass = PhysicsConstants::DEUTERIUM_MASS;
        p.charge = PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 0.3f;
        p.g = 0.6f;
        p.b = 1.0f;
        p.a = 0.9f; // Light blue
        break;
    case Particle::TRITIUM:
        p.mass = PhysicsConstants::TRITIUM_MASS;
        p.charge = PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 0.6f;
        p.g = 0.3f;
        p.b = 1.0f;
        p.a = 0.9f; // Purple
        break;
    case Particle::HELIUM:
        p.mass = PhysicsConstants::HELIUM_MASS;
        p.charge = 2.0f * PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 1.0f;
        p.g = 1.0f;
        p.b = 0.3f;
        p.a = 1.0f; // Yellow
        p.radius = 0.025f;
        break;
    case Particle::NEUTRON:
        p.mass = PhysicsConstants::NEUTRON_MASS;
        p.charge = 0.0f;
        p.r = 0.8f;
        p.g = 0.8f;
        p.b = 0.8f;
        p.a = 0.7f; // Gray
        p.radius = 0.015f;
        break;
    case Particle::ELECTRON:
        p.mass = PhysicsConstants::ELECTRON_MASS;
        p.charge = -PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 1.0f;
        p.g = 0.2f;
        p.b = 0.2f;
        p.a = 0.6f;        // Red
        p.radius = 0.008f; // Even smaller
        break;
    }

    p.kineticEnergy = 0.5f * p.mass * (vx * vx + vy * vy + vz * vz);

    return p;
}

// Convert CPU particle to GPU-uploadable format
inline GPUParticle toGPUParticle(const Particle& p)
{
    GPUParticle gp;
    gp.px = p.x;
    gp.py = p.y;
    gp.pz = p.z;
    gp.radius = p.radius;
    gp.r = p.r;
    gp.g = p.g;
    gp.b = p.b;
    gp.a = p.a;
    return gp;
}

#endif // PARTICLE_H
