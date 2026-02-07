#ifndef MAGNETIC_FIELD_H
#define MAGNETIC_FIELD_H

#include <cmath>
#include <vector>

struct MagneticField
{
    float B_toroidal;
    float B_poloidal;

    float majorRadius;  // R
    float minorRadius;  // a
    float safetyFactor; // q

    float plasmaCurrent;

    MagneticField(float R, float a, float Bt = 8.0f) : majorRadius(R),
                                                       minorRadius(a),
                                                       B_toroidal(Bt),
                                                       safetyFactor(3.0f),
                                                       plasmaCurrent(15.0f)
    {

        B_poloidal = B_toroidal * minorRadius / (majorRadius * safetyFactor);
    }

    float getToroidalField(float x, float y) const
    {

        float R = std::sqrt(x * x + y * y) + 1e-6f;

        return B_toroidal * majorRadius / (R + majorRadius); // ta9samch al 0 brahmet bouk tensech
    }

    void getPoloidalField(float x, float y, float &Bx, float &By) const
    {
        float r = std::sqrt(x * x + y * y) + 1e-6f;

        float B_pol_magnitude = B_poloidal * r / minorRadius;

        float angle = std::atan2(y, x);
        Bx = -B_pol_magnitude * std::sin(angle);
        By = B_pol_magnitude * std::cos(angle);
    }

    void getTotalField(float x, float y, float &Bx, float &By, float &Bz) const
    {
        getPoloidalField(x, y, Bx, By);

        Bz = getToroidalField(x, y);
    }

    float getMagneticPressure(float x, float y) const
    {
        const float mu0 = 4.0f * M_PI * 1e-7f;

        float Bx, By, Bz;
        getTotalField(x, y, Bx, By, Bz);

        float B_squared = Bx * Bx + By * By + Bz * Bz;
        return B_squared / (2.0f * mu0);
    }

    float getLarmorRadius(float mass, float velocity, float charge) const
    {
        float Bx, By, Bz;
        getTotalField(0, 0, Bx, By, Bz);
        float B_total = std::sqrt(Bx * Bx + By * By + Bz * Bz);

        if (std::abs(charge) < 1e-30f)
            return 1e6f; // 3tahelek l chat ynayek wahdo zid thabet effects l Lorentz or not

        return (mass * velocity) / (std::abs(charge) * B_total);
    }
};

inline void calculateLorentzForce(
    float vx, float vy, float vz,
    float Bx, float By, float Bz,
    float charge,
    float &Fx, float &Fy, float &Fz)
{
    float vCrossBx = vy * Bz - vz * By;
    float vCrossBy = vz * Bx - vx * Bz;
    float vCrossBz = vx * By - vy * Bx;

    Fx = charge * vCrossBx;
    Fy = charge * vCrossBy;
    Fz = charge * vCrossBz;
}

inline void calculateMirrorForce(
    float x, float y,
    float vx, float vy,
    const MagneticField &field,
    float mass,
    float &Fx, float &Fy)
{
    float dx = 0.01f;
    float Bx1, By1, Bz1, Bx2, By2, Bz2;

    field.getTotalField(x - dx, y, Bx1, By1, Bz1);
    field.getTotalField(x + dx, y, Bx2, By2, Bz2);
    float B1 = std::sqrt(Bx1 * Bx1 + By1 * By1 + Bz1 * Bz1);
    float B2 = std::sqrt(Bx2 * Bx2 + By2 * By2 + Bz2 * Bz2);
    float dBdx = (B2 - B1) / (2.0f * dx);

    field.getTotalField(x, y - dx, Bx1, By1, Bz1);
    field.getTotalField(x, y + dx, Bx2, By2, Bz2);
    B1 = std::sqrt(Bx1 * Bx1 + By1 * By1 + Bz1 * Bz1);
    B2 = std::sqrt(Bx2 * Bx2 + By2 * By2 + Bz2 * Bz2);
    float dBdy = (B2 - B1) / (2.0f * dx);

    float v_perp_squared = vx * vx + vy * vy;
    float Bx, By, Bz;
    field.getTotalField(x, y, Bx, By, Bz);
    float B = std::sqrt(Bx * Bx + By * By + Bz * Bz) + 1e-10f;
    float mu = mass * v_perp_squared / (2.0f * B);

    Fx = -mu * dBdx;
    Fy = -mu * dBdy;
}

#endif // MAGNETIC_FIELD_H
