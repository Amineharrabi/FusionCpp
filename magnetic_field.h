#ifndef MAGNETIC_FIELD_H
#define MAGNETIC_FIELD_H

#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/**
 * MAGNETIC CONFINEMENT PHYSICS — 3D
 * 
 * Magnetic fields in a Tokamak torus:
 * 1. Toroidal field (Bφ) — wraps around the torus the long way (along the ring)
 * 2. Poloidal field (Bθ) — wraps around the tube cross-section
 * 3. Combined helical field lines confine particles inside the torus tube
 * 
 * Coordinate convention: torus center ring lies in the XZ plane at y=0.
 *   - "toroidal angle" φ: angle around the ring in XZ
 *   - "poloidal angle" θ: angle around the tube cross-section
 */

struct MagneticField {
    // Field strengths (Tesla)
    float B_toroidal;    // Toroidal field strength
    float B_poloidal;    // Poloidal field strength

    // Tokamak geometry parameters
    float majorRadius;   // R
    float minorRadius;   // a
    float safetyFactor;  // q

    // Plasma parameters
    float plasmaCurrent;

    MagneticField(float R, float a, float Bt = 8.0f) :
        majorRadius(R),
        minorRadius(a),
        B_toroidal(Bt),
        safetyFactor(3.0f),
        plasmaCurrent(15.0f)
    {
        B_poloidal = B_toroidal * minorRadius / (majorRadius * safetyFactor);
    }

    /**
     * Get the toroidal field magnitude at a 3D point.
     * The toroidal field falls off as 1/R where R = distance from the torus axis (Y-axis).
     * B_φ = B0 * R0 / R_local
     */
    float getToroidalFieldMagnitude(float x, float y, float z) const {
        float R_local = std::sqrt(x * x + z * z);
        if (R_local < 1e-6f) R_local = 1e-6f;
        return B_toroidal * majorRadius / R_local;
    }

    /**
     * Get the 3D toroidal field direction.
     * The toroidal field circulates around the torus ring axis (Y-axis).
     * At point (x, 0, z), the toroidal direction is (-z, 0, x)/|xz| (tangent to circle)
     */
    void getToroidalFieldDir(float x, float z, float& dx, float& dy, float& dz) const {
        float R = std::sqrt(x * x + z * z);
        if (R < 1e-6f) {
            dx = 0.0f; dy = 0.0f; dz = 1.0f;
            return;
        }
        dx = -z / R;
        dy = 0.0f;
        dz = x / R;
    }

   
    void getPoloidalField3D(float x, float y, float z, float& Bx, float& By, float& Bz) const {
        float Rxz = std::sqrt(x * x + z * z);
        if (Rxz < 1e-6f) Rxz = 1e-6f;
        
        float cx = majorRadius * (x / Rxz);
        float cz = majorRadius * (z / Rxz);
        
        float rx = x - cx;
        float ry = y;  
        float rz = z - cz;
        float rLen = std::sqrt(rx * rx + ry * ry + rz * rz);
        if (rLen < 1e-6f) rLen = 1e-6f;
        
        float rnx = rx / rLen;
        float rny = ry / rLen;
        float rnz = rz / rLen;
        
        float tdx, tdy, tdz;
        getToroidalFieldDir(x, z, tdx, tdy, tdz);
        
      
        float pdx = tdy * rnz - tdz * rny;
        float pdy = tdz * rnx - tdx * rnz;
        float pdz = tdx * rny - tdy * rnx;
        
        float rFrac = rLen / minorRadius;
        if (rFrac > 2.0f) rFrac = 2.0f;
        float B_pol = B_poloidal * rFrac;
        
        Bx = B_pol * pdx;
        By = B_pol * pdy;
        Bz = B_pol * pdz;
    }
    
    void getTotalField(float x, float y, float z, float& Bx, float& By, float& Bz) const {
        // Poloidal field
        getPoloidalField3D(x, y, z, Bx, By, Bz);
        
        // Toroidal field
        float Bt = getToroidalFieldMagnitude(x, y, z);
        float tdx, tdy, tdz;
        getToroidalFieldDir(x, z, tdx, tdy, tdz);
        
        Bx += Bt * tdx;
        By += Bt * tdy;
        Bz += Bt * tdz;
    }

    void getTotalField(float px, float py, float& Bx, float& By, float& Bz) const {
       
        getTotalField(majorRadius + px, py, 0.0f, Bx, By, Bz);
    }

    
    float getMagneticPressure(float x, float y, float z) const {
        const float mu0 = 4.0f * M_PI * 1e-7f;
        float Bx, By, Bz;
        getTotalField(x, y, z, Bx, By, Bz);
        float B_squared = Bx * Bx + By * By + Bz * Bz;
        return B_squared / (2.0f * mu0);
    }

    float getLarmorRadius(float mass, float velocity, float charge) const {
        float Bx, By, Bz;
        getTotalField(majorRadius, 0.0f, 0.0f, Bx, By, Bz);
        float B_total = std::sqrt(Bx * Bx + By * By + Bz * Bz);
        if (std::abs(charge) < 1e-30f) return 1e6f;
        return (mass * velocity) / (std::abs(charge) * B_total);
    }
};


inline void calculateLorentzForce(
    float vx, float vy, float vz,
    float Bx, float By, float Bz,
    float charge,
    float& Fx, float& Fy, float& Fz)
{
    float vCrossBx = vy * Bz - vz * By;
    float vCrossBy = vz * Bx - vx * Bz;
    float vCrossBz = vx * By - vy * Bx;
    
    Fx = charge * vCrossBx;
    Fy = charge * vCrossBy;
    Fz = charge * vCrossBz;
}


inline void calculateMirrorForce3D(
    float x, float y, float z,
    float vx, float vy, float vz,
    const MagneticField& field,
    float mass,
    float& Fx, float& Fy, float& Fz)
{
    const float dx = 0.01f;
    float Bx1, By1, Bz1, Bx2, By2, Bz2;
    
    field.getTotalField(x - dx, y, z, Bx1, By1, Bz1);
    field.getTotalField(x + dx, y, z, Bx2, By2, Bz2);
    float B1 = std::sqrt(Bx1 * Bx1 + By1 * By1 + Bz1 * Bz1);
    float B2 = std::sqrt(Bx2 * Bx2 + By2 * By2 + Bz2 * Bz2);
    float dBdx = (B2 - B1) / (2.0f * dx);
    
    field.getTotalField(x, y - dx, z, Bx1, By1, Bz1);
    field.getTotalField(x, y + dx, z, Bx2, By2, Bz2);
    B1 = std::sqrt(Bx1 * Bx1 + By1 * By1 + Bz1 * Bz1);
    B2 = std::sqrt(Bx2 * Bx2 + By2 * By2 + Bz2 * Bz2);
    float dBdy = (B2 - B1) / (2.0f * dx);
    
    field.getTotalField(x, y, z - dx, Bx1, By1, Bz1);
    field.getTotalField(x, y, z + dx, Bx2, By2, Bz2);
    B1 = std::sqrt(Bx1 * Bx1 + By1 * By1 + Bz1 * Bz1);
    B2 = std::sqrt(Bx2 * Bx2 + By2 * By2 + Bz2 * Bz2);
    float dBdz = (B2 - B1) / (2.0f * dx);
    
    float v_perp_sq = vx * vx + vy * vy + vz * vz;
    float Bx0, By0, Bz0;
    field.getTotalField(x, y, z, Bx0, By0, Bz0);
    float B0 = std::sqrt(Bx0 * Bx0 + By0 * By0 + Bz0 * Bz0) + 1e-10f;
    float mu = mass * v_perp_sq / (2.0f * B0);
    
    Fx = -mu * dBdx;
    Fy = -mu * dBdy;
    Fz = -mu * dBdz;
}

inline void calculateMirrorForce(
    float x, float y,
    float vx, float vy,
    const MagneticField& field,
    float mass,
    float& Fx, float& Fy)
{
    float Fz;
    calculateMirrorForce3D(field.majorRadius + x, y, 0.0f, vx, vy, 0.0f, field, mass, Fx, Fy, Fz);
}

#endif // MAGNETIC_FIELD_H
