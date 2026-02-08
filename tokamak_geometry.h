#ifndef TOKAMAK_GEOMETRY_H
#define TOKAMAK_GEOMETRY_H

#include <vector>
#include <cmath>
#include <glad/glad.h>
#define M_PI 3.141592

/**
 * TOKAMAK GEOMETRY MODULE
 * 
 * This module defines the 2D cross-sectional geometry of a Tokamak reactor.
 * The cross-section shows:
 * - Vacuum vessel (outer containment)
 * - Plasma chamber (where fusion occurs)
 * - First wall (inner boundary)
 * - Divertor regions (for exhaust)
 * - Magnetic coil positions (visual representation)
 */

struct TokamakGeometry {
    // Geometric parameters (in simulation units, mapped to real meters)
    float majorRadius;      // R - distance from center to plasma center (e.g., 2.0m)
    float minorRadius;      // a - plasma cross-section radius (e.g., 0.67m)
    float plasmaElongation; // κ - vertical stretching factor (typical: 1.7-2.0)
    float plasmaTriangularity; // δ - D-shape factor (typical: 0.3-0.5)
    
    // Vessel dimensions
    float vesselThickness;
    float firstWallRadius;
    
    // Rendering data
    std::vector<float> plasmaVertices;
    std::vector<float> vesselVertices;
    std::vector<float> coilVertices;
    std::vector<float> divertorVertices;
    
    GLuint VAO_plasma, VBO_plasma;
    GLuint VAO_vessel, VBO_vessel;
    GLuint VAO_coils, VBO_coils;
    GLuint VAO_divertor, VBO_divertor;
    
    // Initialize with default ITER-like parameters
    TokamakGeometry() :
        majorRadius(0.6f),      // Normalized for screen space
        minorRadius(0.25f),
        plasmaElongation(1.7f),
        plasmaTriangularity(0.33f),
        vesselThickness(0.05f),
        firstWallRadius(0.27f)
    {
        generateGeometry();
        setupBuffers();
    }
    
    void generateGeometry();
    void setupBuffers();
    void render(GLuint shaderProgram);
    void cleanup();
    
    // Check if a point is inside the plasma region
    bool isInsidePlasma(float x, float y) const;
    
    // Get the distance from plasma edge (negative inside, positive outside)
    float distanceFromPlasmaEdge(float x, float y) const;
};

inline float distancePointToSegment(float px, float py, float ax, float ay, float bx, float by)
{
    float abx = bx - ax;
    float aby = by - ay;
    float apx = px - ax;
    float apy = py - ay;

    float abLen2 = abx * abx + aby * aby;
    if (abLen2 < 1e-20f) {
        float dx = px - ax;
        float dy = py - ay;
        return std::sqrt(dx * dx + dy * dy);
    }

    float t = (apx * abx + apy * aby) / abLen2;
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;

    float cx = ax + t * abx;
    float cy = ay + t * aby;
    float dx = px - cx;
    float dy = py - cy;
    return std::sqrt(dx * dx + dy * dy);
}

inline bool pointInPolygonEvenOdd(const std::vector<float>& polyXY, float x, float y)
{
    const size_t n = polyXY.size() / 2;
    if (n < 3) return false;

    bool inside = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        float xi = polyXY[i * 2 + 0];
        float yi = polyXY[i * 2 + 1];
        float xj = polyXY[j * 2 + 0];
        float yj = polyXY[j * 2 + 1];

        bool intersect = ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-20f) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
}

inline float signedDistanceToPolygon(const std::vector<float>& polyXY, float x, float y)
{
    const size_t n = polyXY.size() / 2;
    if (n < 2) return 1e6f;

    float minDist = 1e6f;
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        float ax = polyXY[i * 2 + 0];
        float ay = polyXY[i * 2 + 1];
        float bx = polyXY[j * 2 + 0];
        float by = polyXY[j * 2 + 1];
        float d = distancePointToSegment(x, y, ax, ay, bx, by);
        if (d < minDist) minDist = d;
    }

    bool inside = pointInPolygonEvenOdd(polyXY, x, y);
    return inside ? -minDist : minDist;
}

/**
 * Generate the D-shaped plasma cross-section
 * Uses parametric equations for elongated, triangular plasma
 */
inline void TokamakGeometry::generateGeometry()
{
    plasmaVertices.clear();
    vesselVertices.clear();
    coilVertices.clear();
    divertorVertices.clear();
    
    int segments = 100;
    
    // ==================== PLASMA BOUNDARY ====================
    // Parametric D-shape: combines circular and triangular deformation
    // x(θ) = R + a·cos(θ + δ·sin(θ))
    // y(θ) = κ·a·sin(θ)
    
    for (int i = 0; i <= segments; ++i) {
        float theta = 2.0f * M_PI * i / segments;
        
        // D-shape formula with triangularity
        float r = minorRadius * std::cos(theta + plasmaTriangularity * std::sin(theta));
        float z = plasmaElongation * minorRadius * std::sin(theta);
        
        plasmaVertices.push_back(r);
        plasmaVertices.push_back(z);
    }
    
    // ==================== VESSEL (OUTER WALL) ====================
    float vesselMinorRadius = minorRadius + vesselThickness + 0.1f;
    float vesselElongation = plasmaElongation * 1.1f;
    
    for (int i = 0; i <= segments; ++i) {
        float theta = 2.0f * M_PI * i / segments;
        float r = vesselMinorRadius * std::cos(theta);
        float z = vesselElongation * vesselMinorRadius * std::sin(theta);
        
        vesselVertices.push_back(r);
        vesselVertices.push_back(z);
    }
    
    // ==================== MAGNETIC COILS ====================
    // Toroidal field coils (simplified as rectangles in cross-section)
    int numCoils = 8;
    float coilDistance = vesselMinorRadius + 0.15f;
    float coilWidth = 0.04f;
    float coilHeight = 0.12f;
    
    for (int i = 0; i < numCoils; ++i) {
        float angle = 2.0f * M_PI * i / numCoils;
        float cx = coilDistance * std::cos(angle);
        float cy = coilDistance * std::sin(angle) * vesselElongation;
        
        // Create rectangle for coil
        coilVertices.push_back(cx - coilWidth);
        coilVertices.push_back(cy - coilHeight);
        
        coilVertices.push_back(cx + coilWidth);
        coilVertices.push_back(cy - coilHeight);
        
        coilVertices.push_back(cx + coilWidth);
        coilVertices.push_back(cy + coilHeight);
        
        coilVertices.push_back(cx - coilWidth);
        coilVertices.push_back(cy + coilHeight);
    }
    
    // ==================== DIVERTOR ====================
    // Lower divertor region (simplified)
    float divertorY = -vesselElongation * vesselMinorRadius * 0.9f;
    divertorVertices = {
        -0.15f, divertorY,
        0.15f, divertorY,
        0.12f, divertorY - 0.08f,
        -0.12f, divertorY - 0.08f
    };
}

inline void TokamakGeometry::setupBuffers()
{
    // Plasma boundary
    glGenVertexArrays(1, &VAO_plasma);
    glGenBuffers(1, &VBO_plasma);
    glBindVertexArray(VAO_plasma);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_plasma);
    glBufferData(GL_ARRAY_BUFFER, plasmaVertices.size() * sizeof(float), 
                 plasmaVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Vessel
    glGenVertexArrays(1, &VAO_vessel);
    glGenBuffers(1, &VBO_vessel);
    glBindVertexArray(VAO_vessel);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_vessel);
    glBufferData(GL_ARRAY_BUFFER, vesselVertices.size() * sizeof(float), 
                 vesselVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Coils
    glGenVertexArrays(1, &VAO_coils);
    glGenBuffers(1, &VBO_coils);
    glBindVertexArray(VAO_coils);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_coils);
    glBufferData(GL_ARRAY_BUFFER, coilVertices.size() * sizeof(float), 
                 coilVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Divertor
    glGenVertexArrays(1, &VAO_divertor);
    glGenBuffers(1, &VBO_divertor);
    glBindVertexArray(VAO_divertor);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_divertor);
    glBufferData(GL_ARRAY_BUFFER, divertorVertices.size() * sizeof(float), 
                 divertorVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    glBindVertexArray(0);
}

inline void TokamakGeometry::render(GLuint shaderProgram)
{
    GLint colorLoc = glGetUniformLocation(shaderProgram, "particleColor");
    
    // Draw vessel (dark gray)
    glBindVertexArray(VAO_vessel);
    glUniform4f(colorLoc, 0.2f, 0.2f, 0.25f, 1.0f);
    glDrawArrays(GL_LINE_LOOP, 0, vesselVertices.size() / 2);
    
    // Draw plasma boundary (cyan glow)
    glBindVertexArray(VAO_plasma);
    glUniform4f(colorLoc, 0.0f, 0.8f, 1.0f, 0.6f);
    glDrawArrays(GL_LINE_LOOP, 0, plasmaVertices.size() / 2);
    
    // Draw magnetic coils (orange/copper)
    glBindVertexArray(VAO_coils);
    glUniform4f(colorLoc, 0.9f, 0.5f, 0.1f, 1.0f);
    for (size_t i = 0; i < coilVertices.size() / 8; ++i) {
        glDrawArrays(GL_LINE_LOOP, i * 4, 4);
    }
    
    // Draw divertor (red)
    glBindVertexArray(VAO_divertor);
    glUniform4f(colorLoc, 0.8f, 0.2f, 0.2f, 1.0f);
    glDrawArrays(GL_LINE_LOOP, 0, divertorVertices.size() / 2);
    
    glBindVertexArray(0);
}

inline bool TokamakGeometry::isInsidePlasma(float x, float y) const
{
    return signedDistanceToPolygon(plasmaVertices, x, y) <= 0.0f;
}

inline float TokamakGeometry::distanceFromPlasmaEdge(float x, float y) const
{
    return signedDistanceToPolygon(plasmaVertices, x, y);
}

inline void TokamakGeometry::cleanup()
{
    glDeleteVertexArrays(1, &VAO_plasma);
    glDeleteBuffers(1, &VBO_plasma);
    glDeleteVertexArrays(1, &VAO_vessel);
    glDeleteBuffers(1, &VBO_vessel);
    glDeleteVertexArrays(1, &VAO_coils);
    glDeleteBuffers(1, &VBO_coils);
    glDeleteVertexArrays(1, &VAO_divertor);
    glDeleteBuffers(1, &VBO_divertor);
}

#endif // TOKAMAK_GEOMETRY_H
