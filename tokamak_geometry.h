#ifndef TOKAMAK_GEOMETRY_H
#define TOKAMAK_GEOMETRY_H

#include <vector>
#include <cmath>
#include <glad/glad.h>
#define M_PI 3.141592

struct TokamakGeometry
{
    float majorRadius;         // R
    float minorRadius;         // a
    float plasmaElongation;    // κ
    float plasmaTriangularity; // δ

    float vesselThickness;
    float firstWallRadius;

    std::vector<float> plasmaVertices;
    std::vector<float> vesselVertices;
    std::vector<float> coilVertices;
    std::vector<float> divertorVertices;

    GLuint VAO_plasma, VBO_plasma;
    GLuint VAO_vessel, VBO_vessel;
    GLuint VAO_coils, VBO_coils;
    GLuint VAO_divertor, VBO_divertor;

    // thabet ml ITER papers
    TokamakGeometry() : majorRadius(0.6f),
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

    bool isInsidePlasma(float x, float y) const;

    float distanceFromPlasmaEdge(float x, float y) const;
};

inline void TokamakGeometry::generateGeometry()
{
    plasmaVertices.clear();
    vesselVertices.clear();
    coilVertices.clear();
    divertorVertices.clear();

    int segments = 100;

    for (int i = 0; i <= segments; ++i)
    {
        float theta = 2.0f * M_PI * i / segments;

        float r = minorRadius * std::cos(theta + plasmaTriangularity * std::sin(theta));
        float z = plasmaElongation * minorRadius * std::sin(theta);

        plasmaVertices.push_back(r);
        plasmaVertices.push_back(z);
    }

    float vesselMinorRadius = minorRadius + vesselThickness + 0.1f;
    float vesselElongation = plasmaElongation * 1.1f;

    for (int i = 0; i <= segments; ++i)
    {
        float theta = 2.0f * M_PI * i / segments;
        float r = vesselMinorRadius * std::cos(theta);
        float z = vesselElongation * vesselMinorRadius * std::sin(theta);

        vesselVertices.push_back(r);
        vesselVertices.push_back(z);
    }

    int numCoils = 8;
    float coilDistance = vesselMinorRadius + 0.15f;
    float coilWidth = 0.04f;
    float coilHeight = 0.12f;

    for (int i = 0; i < numCoils; ++i)
    {
        float angle = 2.0f * M_PI * i / numCoils;
        float cx = coilDistance * std::cos(angle);
        float cy = coilDistance * std::sin(angle) * vesselElongation;

        coilVertices.push_back(cx - coilWidth);
        coilVertices.push_back(cy - coilHeight);

        coilVertices.push_back(cx + coilWidth);
        coilVertices.push_back(cy - coilHeight);

        coilVertices.push_back(cx + coilWidth);
        coilVertices.push_back(cy + coilHeight);

        coilVertices.push_back(cx - coilWidth);
        coilVertices.push_back(cy + coilHeight);
    }

    float divertorY = -vesselElongation * vesselMinorRadius * 0.9f;
    divertorVertices = {
        -0.15f, divertorY,
        0.15f, divertorY,
        0.12f, divertorY - 0.08f,
        -0.12f, divertorY - 0.08f};
}

inline void TokamakGeometry::setupBuffers()
{
    glGenVertexArrays(1, &VAO_plasma);
    glGenBuffers(1, &VBO_plasma);
    glBindVertexArray(VAO_plasma);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_plasma);
    glBufferData(GL_ARRAY_BUFFER, plasmaVertices.size() * sizeof(float),
                 plasmaVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    glGenVertexArrays(1, &VAO_vessel);
    glGenBuffers(1, &VBO_vessel);
    glBindVertexArray(VAO_vessel);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_vessel);
    glBufferData(GL_ARRAY_BUFFER, vesselVertices.size() * sizeof(float),
                 vesselVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    glGenVertexArrays(1, &VAO_coils);
    glGenBuffers(1, &VBO_coils);
    glBindVertexArray(VAO_coils);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_coils);
    glBufferData(GL_ARRAY_BUFFER, coilVertices.size() * sizeof(float),
                 coilVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    glGenVertexArrays(1, &VAO_divertor);
    glGenBuffers(1, &VBO_divertor);
    glBindVertexArray(VAO_divertor);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_divertor);
    glBufferData(GL_ARRAY_BUFFER, divertorVertices.size() * sizeof(float),
                 divertorVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}

inline void TokamakGeometry::render(GLuint shaderProgram)
{
    GLint colorLoc = glGetUniformLocation(shaderProgram, "particleColor");

    glBindVertexArray(VAO_vessel);
    glUniform4f(colorLoc, 0.2f, 0.2f, 0.25f, 1.0f);
    glDrawArrays(GL_LINE_LOOP, 0, vesselVertices.size() / 2);

    glBindVertexArray(VAO_plasma);
    glUniform4f(colorLoc, 0.0f, 0.8f, 1.0f, 0.6f);
    glDrawArrays(GL_LINE_LOOP, 0, plasmaVertices.size() / 2);

    glBindVertexArray(VAO_coils);
    glUniform4f(colorLoc, 0.9f, 0.5f, 0.1f, 1.0f);
    for (size_t i = 0; i < coilVertices.size() / 8; ++i)
    {
        glDrawArrays(GL_LINE_LOOP, i * 4, 4);
    }

    glBindVertexArray(VAO_divertor);
    glUniform4f(colorLoc, 0.8f, 0.2f, 0.2f, 1.0f);
    glDrawArrays(GL_LINE_LOOP, 0, divertorVertices.size() / 2);

    glBindVertexArray(0);
}

inline bool TokamakGeometry::isInsidePlasma(float x, float y) const
{

    float normalizedX = x / minorRadius;
    float normalizedY = y / (plasmaElongation * minorRadius);

    float distance = normalizedX * normalizedX + normalizedY * normalizedY;

    return distance <= 1.0f;
}

inline float TokamakGeometry::distanceFromPlasmaEdge(float x, float y) const
{
    float normalizedX = x / minorRadius;
    float normalizedY = y / (plasmaElongation * minorRadius);
    float distance = std::sqrt(normalizedX * normalizedX + normalizedY * normalizedY);

    return distance - 1.0f;
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
