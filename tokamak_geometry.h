#ifndef TOKAMAK_GEOMETRY_H
#define TOKAMAK_GEOMETRY_H

#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/**
 * TOKAMAK GEOMETRY MODULE
 * 
 * This module defines both 2D cross-sectional and 3D toroidal geometry
 * of a Tokamak reactor.
 * 
 * 3D Torus: the plasma volume is a torus with major radius R and minor radius r.
 * A point (x, y, z) is inside the torus when:
 *   (sqrt(x² + z²) - R)² + y² < r²
 * 
 * The torus center ring lies in the XZ plane at y=0.
 */

struct TokamakGeometry {
    // 3D Torus parameters (simulation units)
    float torusMajorR;       // R - major radius of the torus
    float torusMinorR;       // r - minor radius (tube thickness)
    float torusOpacity;      // Rendering opacity for the shell

    // 2D cross-section parameters (kept for physics compatibility)
    float majorRadius;       // R (same as torusMajorR)
    float minorRadius;       // a (same as torusMinorR)
    float plasmaElongation;
    float plasmaTriangularity;
    float vesselThickness;
    float firstWallRadius;
    
    // 2D rendering data (kept for reference / debug overlay)
    std::vector<float> plasmaVertices;
    std::vector<float> vesselVertices;

    // Initialize with default ITER-like parameters
    TokamakGeometry() :
        torusMajorR(1.2f),
        torusMinorR(0.4f),
        torusOpacity(0.15f),
        majorRadius(1.2f),
        minorRadius(0.4f),
        plasmaElongation(1.7f),
        plasmaTriangularity(0.33f),
        vesselThickness(0.05f),
        firstWallRadius(0.42f)
    {
        generateCrossSection();
    }

    // ==================== 3D TORUS SDF ====================
    
    /**
     * Signed Distance Function for a torus centered at origin in the XZ plane.
     * Negative = inside, Positive = outside.
     * SDF(p) = length(vec2(length(p.xz) - R, p.y)) - r
     */
    float torusSDF(float x, float y, float z) const {
        float dxz = std::sqrt(x * x + z * z) - torusMajorR;
        return std::sqrt(dxz * dxz + y * y) - torusMinorR;
    }
    
    /**
     * Check if a 3D point is inside the plasma torus volume
     */
    bool isInsidePlasma3D(float x, float y, float z) const {
        return torusSDF(x, y, z) <= 0.0f;
    }
    
    /**
     * Get signed distance from the torus surface (negative inside)
     */
    float distanceFromPlasmaEdge3D(float x, float y, float z) const {
        return torusSDF(x, y, z);
    }

    /**
     * Get the gradient of the torus SDF (outward normal direction)
     */
    void torusNormal(float x, float y, float z, float& nx, float& ny, float& nz) const {
        const float eps = 0.001f;
        float d = torusSDF(x, y, z);
        nx = torusSDF(x + eps, y, z) - d;
        ny = torusSDF(x, y + eps, z) - d;
        nz = torusSDF(x, y, z + eps) - d;
        float len = std::sqrt(nx * nx + ny * ny + nz * nz) + 1e-10f;
        nx /= len;
        ny /= len;
        nz /= len;
    }
    
    /**
     * Project a point onto the nearest point on the torus centerline ring.
     * The centerline is a circle of radius R in the XZ plane at y=0.
     * Returns the point on the centerline closest to (x, y, z).
     */
    void projectToCenterline(float x, float y, float z,
                             float& cx, float& cy, float& cz) const {
        float rxz = std::sqrt(x * x + z * z);
        if (rxz < 1e-8f) {
            // On the torus axis — pick an arbitrary point on the centerline
            cx = torusMajorR;
            cy = 0.0f;
            cz = 0.0f;
        } else {
            cx = torusMajorR * (x / rxz);
            cy = 0.0f;
            cz = torusMajorR * (z / rxz);
        }
    }

    // ==================== 2D cross-section (reference) ====================

    /**
     * Generate 2D D-shaped plasma cross-section (for ImGui overlay / debug)
     */
    void generateCrossSection() {
        plasmaVertices.clear();
        vesselVertices.clear();
        int segments = 100;

        // Plasma boundary (D-shape)
        for (int i = 0; i <= segments; ++i) {
            float theta = 2.0f * M_PI * i / segments;
            float r = minorRadius * std::cos(theta + plasmaTriangularity * std::sin(theta));
            float zz = plasmaElongation * minorRadius * std::sin(theta);
            plasmaVertices.push_back(r);
            plasmaVertices.push_back(zz);
        }

        // Vessel boundary
        float vesselR = minorRadius + vesselThickness + 0.1f;
        float vesselE = plasmaElongation * 1.1f;
        for (int i = 0; i <= segments; ++i) {
            float theta = 2.0f * M_PI * i / segments;
            float r = vesselR * std::cos(theta);
            float zz = vesselE * vesselR * std::sin(theta);
            vesselVertices.push_back(r);
            vesselVertices.push_back(zz);
        }
    }

    // 2D containment check (legacy — uses cross-section polygon)
    bool isInsidePlasma(float x, float y) const {
        // Simple elliptical check for the 2D cross-section
        float dxz = x; // In 2D cross-section, x represents radial offset from center
        float normalizedR = dxz / minorRadius;
        float normalizedZ = y / (plasmaElongation * minorRadius);
        return (normalizedR * normalizedR + normalizedZ * normalizedZ) <= 1.0f;
    }

    float distanceFromPlasmaEdge(float x, float y) const {
        float normalizedR = x / minorRadius;
        float normalizedZ = y / (plasmaElongation * minorRadius);
        float ellipseVal = normalizedR * normalizedR + normalizedZ * normalizedZ;
        // Approximate signed distance
        float scale = (minorRadius + plasmaElongation * minorRadius) * 0.5f;
        return (std::sqrt(ellipseVal) - 1.0f) * scale;
    }

    void cleanup() {
        // No GPU resources in 3D mode (compute shader handles rendering)
    }
};

#endif // TOKAMAK_GEOMETRY_H
