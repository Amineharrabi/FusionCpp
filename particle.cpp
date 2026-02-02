#include "particle.h"
#include <cmath>

std::vector<float> generateCircleVertices(float cx, float cy, float radius, int segments)
{
    std::vector<float> vertices;
    float angleStep = 2.0f * 3.1415926f / segments;
    // Center vertex for triangle fan
    vertices.push_back(cx);
    vertices.push_back(cy);
    for (int i = 0; i <= segments; ++i)
    {
        float angle = i * angleStep;
        float x = cx + radius * cosf(angle);
        float y = cy + radius * sinf(angle);
        vertices.push_back(x);
        vertices.push_back(y);
    }
    return vertices;
}
