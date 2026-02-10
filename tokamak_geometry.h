#ifndef TOKAMAK_GEOMETRY_H
#define TOKAMAK_GEOMETRY_H

#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif



struct TokamakGeometry {
    
    float torusMajorR;       
    float torusMinorR;       
    float torusOpacity;     

    float majorRadius;       
    float minorRadius;       
    float plasmaElongation;
    float plasmaTriangularity;
    float vesselThickness;
    float firstWallRadius;
    
    std::vector<float> plasmaVertices;
    std::vector<float> vesselVertices;

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

   
    float torusSDF(float x, float y, float z) const {
        float dxz = std::sqrt(x * x + z * z) - torusMajorR;
        return std::sqrt(dxz * dxz + y * y) - torusMinorR;
    }
    
    
    bool isInsidePlasma3D(float x, float y, float z) const {
        return torusSDF(x, y, z) <= 0.0f;
    }
    
    
    float distanceFromPlasmaEdge3D(float x, float y, float z) const {
        return torusSDF(x, y, z);
    }

   
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
    
   
    void projectToCenterline(float x, float y, float z,
                             float& cx, float& cy, float& cz) const {
        float rxz = std::sqrt(x * x + z * z);
        if (rxz < 1e-8f) {
            cx = torusMajorR;
            cy = 0.0f;
            cz = 0.0f;
        } else {
            cx = torusMajorR * (x / rxz);
            cy = 0.0f;
            cz = torusMajorR * (z / rxz);
        }
    }

    
    void generateCrossSection() {
        plasmaVertices.clear();
        vesselVertices.clear();
        int segments = 100;

        for (int i = 0; i <= segments; ++i) {
            float theta = 2.0f * M_PI * i / segments;
            float r = minorRadius * std::cos(theta + plasmaTriangularity * std::sin(theta));
            float zz = plasmaElongation * minorRadius * std::sin(theta);
            plasmaVertices.push_back(r);
            plasmaVertices.push_back(zz);
        }

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

    bool isInsidePlasma(float x, float y) const {
        float dxz = x; 
        float normalizedR = dxz / minorRadius;
        float normalizedZ = y / (plasmaElongation * minorRadius);
        return (normalizedR * normalizedR + normalizedZ * normalizedZ) <= 1.0f;
    }

    float distanceFromPlasmaEdge(float x, float y) const {
        float normalizedR = x / minorRadius;
        float normalizedZ = y / (plasmaElongation * minorRadius);
        float ellipseVal = normalizedR * normalizedR + normalizedZ * normalizedZ;
        float scale = (minorRadius + plasmaElongation * minorRadius) * 0.5f;
        return (std::sqrt(ellipseVal) - 1.0f) * scale;
    }

    void cleanup() {
    }
};

#endif // TOKAMAK_GEOMETRY_H
