#ifndef MAGNETIC_FIELD_H
#define MAGNETIC_FIELD_H

#include <cmath>
#include <vector>

/**
 * MAGNETIC CONFINEMENT PHYSICS
 * 
 * This module calculates magnetic fields in a Tokamak:
 * 1. Toroidal field (Bφ) - main confining field around the torus
 * 2. Poloidal field (Bθ) - field around the plasma cross-section
 * 3. Total field and its effects on charged particles
 */

struct MagneticField {
    // Field strengths (Tesla)
    float B_toroidal;    // Toroidal field strength (typical: 5-15 T)
    float B_poloidal;    // Poloidal field strength (typical: 0.5-2 T)
    
    // Tokamak geometry parameters
    float majorRadius;   // R
    float minorRadius;   // a
    float safetyFactor;  // q - magnetic safety factor (typical: 1-4)
    
    // Plasma parameters
    float plasmaCurrent; // Ip (Mega-Amperes in real tokamaks)
    
    MagneticField(float R, float a, float Bt = 8.0f) :
        majorRadius(R),
        minorRadius(a),
        B_toroidal(Bt),
        safetyFactor(3.0f),
        plasmaCurrent(15.0f)
    {
        // Poloidal field calculated from plasma current and safety factor
        // B_θ ≈ μ₀·Ip / (2π·r)
        B_poloidal = B_toroidal * minorRadius / (majorRadius * safetyFactor);
    }
    
    /**
     * Calculate toroidal magnetic field at position (x, y)
     * B_φ = B₀·R₀/R where R is distance from major axis
     * 
     * In 2D cross-section, this field points into/out of screen
     * but affects particle motion through Lorentz force
     */
    float getToroidalField(float x, float y) const {
        // Distance from center (approximation for 2D)
        float R = std::sqrt(x*x + y*y) + 1e-6f; // Avoid division by zero
        
        // Toroidal field falls off as 1/R
        return B_toroidal * majorRadius / (R + majorRadius);
    }
    
    /**
     * Calculate poloidal magnetic field components at (x, y)
     * The poloidal field circulates around the plasma cross-section
     * 
     * B_poloidal = (Bx, By) in the cross-sectional plane
     */
    void getPoloidalField(float x, float y, float& Bx, float& By) const {
        // Distance from plasma center
        float r = std::sqrt(x*x + y*y) + 1e-6f;
        
        // Poloidal field strength increases linearly with radius (simplified)
        float B_pol_magnitude = B_poloidal * r / minorRadius;
        
        // Field direction: perpendicular to radius (tangential)
        // This creates circulation around the plasma
        float angle = std::atan2(y, x);
        Bx = -B_pol_magnitude * std::sin(angle); // Perpendicular component
        By = B_pol_magnitude * std::cos(angle);
    }
    
    /**
     * Calculate total magnetic field at position (x, y)
     * Returns field components in 2D plane
     */
    void getTotalField(float x, float y, float& Bx, float& By, float& Bz) const {
        // Poloidal field (in-plane)
        getPoloidalField(x, y, Bx, By);
        
        // Toroidal field (out-of-plane, but affects dynamics)
        Bz = getToroidalField(x, y);
    }
    
    /**
     * Calculate magnetic pressure at position
     * P_mag = B²/(2μ₀)
     */
    float getMagneticPressure(float x, float y) const {
        const float mu0 = 4.0f * M_PI * 1e-7f; // Permeability of free space
        
        float Bx, By, Bz;
        getTotalField(x, y, Bx, By, Bz);
        
        float B_squared = Bx*Bx + By*By + Bz*Bz;
        return B_squared / (2.0f * mu0);
    }
    
    /**
     * Calculate Larmor radius (gyroradius) for a particle
     * r_L = m·v_perp / (q·B)
     * 
     * This is the radius of circular motion of charged particle in B-field
     */
    float getLarmorRadius(float mass, float velocity, float charge) const {
        float Bx, By, Bz;
        getTotalField(0, 0, Bx, By, Bz); // Use center field as reference
        float B_total = std::sqrt(Bx*Bx + By*By + Bz*Bz);
        
        if (std::abs(charge) < 1e-30f) return 1e6f; // Neutral particles not affected
        
        return (mass * velocity) / (std::abs(charge) * B_total);
    }
};

/**
 * Calculate Lorentz force on a charged particle
 * F = q(E + v × B)
 * 
 * In Tokamak, E-field is usually small, dominated by v × B force
 */
inline void calculateLorentzForce(
    float vx, float vy, float vz,  // Particle velocity
    float Bx, float By, float Bz,  // Magnetic field
    float charge,                   // Particle charge
    float& Fx, float& Fy, float& Fz) // Output force
{
    // Cross product: v × B
    float vCrossBx = vy * Bz - vz * By;
    float vCrossBy = vz * Bx - vx * Bz;
    float vCrossBz = vx * By - vy * Bx;
    
    // F = q(v × B)
    Fx = charge * vCrossBx;
    Fy = charge * vCrossBy;
    Fz = charge * vCrossBz;
}

/**
 * Magnetic mirror force for particle confinement
 * Particles with velocity component parallel to field lines
 * experience a force when field strength changes
 * 
 * F_mirror = -μ·∇B where μ is magnetic moment
 */
inline void calculateMirrorForce(
    float x, float y,
    float vx, float vy,
    const MagneticField& field,
    float mass,
    float& Fx, float& Fy)
{
    // Calculate magnetic field gradient (simplified)
    float dx = 0.01f;
    float Bx1, By1, Bz1, Bx2, By2, Bz2;
    
    field.getTotalField(x - dx, y, Bx1, By1, Bz1);
    field.getTotalField(x + dx, y, Bx2, By2, Bz2);
    float B1 = std::sqrt(Bx1*Bx1 + By1*By1 + Bz1*Bz1);
    float B2 = std::sqrt(Bx2*Bx2 + By2*By2 + Bz2*Bz2);
    float dBdx = (B2 - B1) / (2.0f * dx);
    
    field.getTotalField(x, y - dx, Bx1, By1, Bz1);
    field.getTotalField(x, y + dx, Bx2, By2, Bz2);
    B1 = std::sqrt(Bx1*Bx1 + By1*By1 + Bz1*Bz1);
    B2 = std::sqrt(Bx2*Bx2 + By2*By2 + Bz2*Bz2);
    float dBdy = (B2 - B1) / (2.0f * dx);
    
    // Magnetic moment: μ = m·v_perp² / (2B)
    float v_perp_squared = vx*vx + vy*vy;
    float Bx, By, Bz;
    field.getTotalField(x, y, Bx, By, Bz);
    float B = std::sqrt(Bx*Bx + By*By + Bz*Bz) + 1e-10f;
    float mu = mass * v_perp_squared / (2.0f * B);
    
    // Mirror force: F = -μ·∇B
    Fx = -mu * dBdx;
    Fy = -mu * dBdy;
}

#endif // MAGNETIC_FIELD_H
