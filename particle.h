#pragma once
#include <vector>

struct Particle
{
    float x, y;       // Position
    float radius;     // Radius
    float r, g, b, a; // Color (RGBA)
};

// Generates (x, y) pairs for a circle centered at (cx, cy) with given radius and segment count
std::vector<float> generateCircleVertices(float cx, float cy, float radius, int segments = 32);