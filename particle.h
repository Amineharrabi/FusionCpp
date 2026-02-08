#ifndef PARTICLE_H
#define PARTICLE_H

#define M_PI 3.141592
#include <vector>
#include <cmath>

struct Particle
{
    // Position (in meters, simulation space)
    float x, y;

    // Velocity (m/s)
    float vx, vy;

    // Physical properties
    float mass;   // kg (e.g., 1.673e-27 for proton, 3.344e-27 for deuterium)
    float charge; // Coulombs (e.g., 1.602e-19 for proton)
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

// Physical constants
namespace PhysicsConstants
{
    constexpr float ELECTRON_MASS = 9.109e-31f;  // kg
    constexpr float PROTON_MASS = 1.673e-27f;    // kg
    constexpr float DEUTERIUM_MASS = 3.344e-27f; // kg (2 × proton mass approx)
    constexpr float TRITIUM_MASS = 5.008e-27f;   // kg (3 × proton mass approx)
    constexpr float HELIUM_MASS = 6.646e-27f;    // kg (4 × proton mass approx)
    constexpr float NEUTRON_MASS = 1.675e-27f;   // kg

    constexpr float ELEMENTARY_CHARGE = 1.602e-19f;   // Coulombs
    constexpr float VACUUM_PERMITTIVITY = 8.854e-12f; // F/m (ε₀)
    constexpr float COULOMB_CONSTANT = 8.988e9f;      // N⋅m²/C² (k = 1/(4πε₀))
    constexpr float BOLTZMANN_CONSTANT = 1.381e-23f;  // J/K

    // Fusion cross-section parameters (simplified)
    constexpr float FUSION_THRESHOLD_ENERGY = 1.0e-14f; // ~60 keV in Joules
    constexpr float FUSION_CROSS_SECTION = 1.0e-28f;    // m² (simplified)
}

// Generate circle vertices for rendering particles
inline std::vector<float> generateCircleVertices(float cx, float cy, float radius, int segments)
{
    std::vector<float> vertices;
    vertices.push_back(cx);
    vertices.push_back(cy);

    for (int i = 0; i <= segments; ++i)
    {
        float angle = 2.0f * M_PI * i / segments;
        vertices.push_back(cx + radius * std::cos(angle));
        vertices.push_back(cy + radius * std::sin(angle));
    }

    return vertices;
}

// Create particle with default properties based on type
inline Particle createParticle(Particle::Type type, float x, float y, float vx, float vy)
{
    Particle p;
    p.x = x;
    p.y = y;
    p.vx = vx;
    p.vy = vy;
    p.type = type;
    p.active = true;
    p.radius = 0.006f; // Very small visual radius

    switch (type)
    {
    case Particle::DEUTERIUM:
        p.mass = PhysicsConstants::DEUTERIUM_MASS;
        p.charge = PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 0.3f;
        p.g = 0.6f;
        p.b = 1.0f;
        p.a = 0.8f; // Light blue
        break;
    case Particle::TRITIUM:
        p.mass = PhysicsConstants::TRITIUM_MASS;
        p.charge = PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 0.6f;
        p.g = 0.3f;
        p.b = 1.0f;
        p.a = 0.8f; // Purple
        break;
    case Particle::HELIUM:
        p.mass = PhysicsConstants::HELIUM_MASS;
        p.charge = 2.0f * PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 1.0f;
        p.g = 1.0f;
        p.b = 0.3f;
        p.a = 0.9f; // Yellow
        break;
    case Particle::NEUTRON:
        p.mass = PhysicsConstants::NEUTRON_MASS;
        p.charge = 0.0f;
        p.r = 0.8f;
        p.g = 0.8f;
        p.b = 0.8f;
        p.a = 0.7f; // Gray
        break;
    case Particle::ELECTRON:
        p.mass = PhysicsConstants::ELECTRON_MASS;
        p.charge = -PhysicsConstants::ELEMENTARY_CHARGE;
        p.r = 1.0f;
        p.g = 0.2f;
        p.b = 0.2f;
        p.a = 0.6f;        // Red
        p.radius = 0.003f; // Even smaller
        break;
    }

    p.kineticEnergy = 0.5f * p.mass * (vx * vx + vy * vy);

    return p;
}

#endif // PARTICLE_H
