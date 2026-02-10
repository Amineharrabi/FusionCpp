#ifndef RAY_TRACING_H
#define RAY_TRACING_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "particle.h"


struct SimulationUBO {
    glm::mat4 invViewProj;    
    glm::vec4 cameraPos;      
    glm::vec4 torusParams;    
    glm::ivec4 counts;        
};

class GPURayTracer {
public:
    GLuint computeProgram = 0;
    GLuint outputTexture = 0;

    GLuint blitProgram = 0;
    GLuint quadVAO = 0;
    GLuint quadVBO = 0;

    GLuint simulationUBO = 0;
    GLuint particleSSBO = 0;
    GLuint flashSSBO = 0;

    int width = 1200;
    int height = 800;

    static const int MAX_PARTICLES = 20000;
    static const int MAX_FLASHES = 64;

    bool initialize(int w, int h) {
        width = w;
        height = h;

        if (!initComputeShader()) return false;
        if (!initBlitShader()) return false;
        createFullscreenQuad();
        createOutputTexture();
        createBuffers();

        std::cout << "GPU Ray Tracer initialized (" << width << "x" << height << ")" << std::endl;
        return true;
    }

    void render(const glm::mat4& invViewProj,
                const glm::vec3& cameraPos,
                float torusMajorR, float torusMinorR, float torusOpacity,
                float time,
                const std::vector<GPUParticle>& gpuParticles,
                const std::vector<FusionFlash>& fusionFlashes,
                int numParticles)
    {
        SimulationUBO ubo;
        ubo.invViewProj = invViewProj;
        ubo.cameraPos = glm::vec4(cameraPos, 0.0f);
        ubo.torusParams = glm::vec4(torusMajorR, torusMinorR, torusOpacity, time);
        ubo.counts = glm::ivec4(
            numParticles,
            (int)fusionFlashes.size(),
            width,
            height
        );

        glBindBuffer(GL_UNIFORM_BUFFER, simulationUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(SimulationUBO), &ubo);

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, particleSSBO);
        if (!gpuParticles.empty()) {
            size_t dataSize = gpuParticles.size() * sizeof(GPUParticle);
            if (dataSize > MAX_PARTICLES * sizeof(GPUParticle))
                dataSize = MAX_PARTICLES * sizeof(GPUParticle);
            glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, dataSize, gpuParticles.data());
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, flashSSBO);
        if (!fusionFlashes.empty()) {
            size_t dataSize = fusionFlashes.size() * sizeof(FusionFlash);
            if (dataSize > MAX_FLASHES * sizeof(FusionFlash))
                dataSize = MAX_FLASHES * sizeof(FusionFlash);
            glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, dataSize, fusionFlashes.data());
        }

        glUseProgram(computeProgram);

        glBindBufferBase(GL_UNIFORM_BUFFER, 0, simulationUBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, particleSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, flashSSBO);
        glBindImageTexture(0, outputTexture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);

        int groupsX = (width + 15) / 16;
        int groupsY = (height + 15) / 16;
        glDispatchCompute(groupsX, groupsY, 1);

        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        blitToScreen();
    }

    void blitToScreen() {
        glUseProgram(blitProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, outputTexture);
        GLint texLoc = glGetUniformLocation(blitProgram, "screenTexture");
        glUniform1i(texLoc, 0);

        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    }

    void resize(int w, int h) {
        width = w;
        height = h;
        if (outputTexture) glDeleteTextures(1, &outputTexture);
        createOutputTexture();
    }

    void cleanup() {
        if (computeProgram) glDeleteProgram(computeProgram);
        if (blitProgram) glDeleteProgram(blitProgram);
        if (quadVAO) glDeleteVertexArrays(1, &quadVAO);
        if (quadVBO) glDeleteBuffers(1, &quadVBO);
        if (outputTexture) glDeleteTextures(1, &outputTexture);
        if (simulationUBO) glDeleteBuffers(1, &simulationUBO);
        if (particleSSBO) glDeleteBuffers(1, &particleSSBO);
        if (flashSSBO) glDeleteBuffers(1, &flashSSBO);
    }

private:
    std::string loadFile(const char* path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "ERROR: Cannot open file: " << path << std::endl;
            return "";
        }
        std::stringstream buf;
        buf << file.rdbuf();
        return buf.str();
    }

    GLuint compileShader(GLenum type, const char* src, const char* label) {
        GLuint shader = glCreateShader(type);
        glShaderSource(shader, 1, &src, nullptr);
        glCompileShader(shader);

        GLint success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char log[2048];
            glGetShaderInfoLog(shader, 2048, nullptr, log);
            std::cerr << "Shader compile error (" << label << "):\n" << log << std::endl;
            return 0;
        }
        return shader;
    }

    bool initComputeShader() {
        std::string compSrc = loadFile("tokamak_raytrace.comp");
        if (compSrc.empty()) {
            std::cerr << "Failed to load compute shader" << std::endl;
            return false;
        }

        GLuint compShader = compileShader(GL_COMPUTE_SHADER, compSrc.c_str(), "compute");
        if (!compShader) return false;

        computeProgram = glCreateProgram();
        glAttachShader(computeProgram, compShader);
        glLinkProgram(computeProgram);

        GLint success;
        glGetProgramiv(computeProgram, GL_LINK_STATUS, &success);
        if (!success) {
            char log[2048];
            glGetProgramInfoLog(computeProgram, 2048, nullptr, log);
            std::cerr << "Compute program link error:\n" << log << std::endl;
            glDeleteShader(compShader);
            return false;
        }

        glDeleteShader(compShader);
        std::cout << "Compute shader compiled and linked successfully" << std::endl;
        return true;
    }

    bool initBlitShader() {
        std::string vertSrc = loadFile("particle.vert");
        std::string fragSrc = loadFile("particle.frag");
        if (vertSrc.empty() || fragSrc.empty()) {
            std::cerr << "Failed to load blit shaders" << std::endl;
            return false;
        }

        GLuint vertShader = compileShader(GL_VERTEX_SHADER, vertSrc.c_str(), "blit_vert");
        GLuint fragShader = compileShader(GL_FRAGMENT_SHADER, fragSrc.c_str(), "blit_frag");
        if (!vertShader || !fragShader) return false;

        blitProgram = glCreateProgram();
        glAttachShader(blitProgram, vertShader);
        glAttachShader(blitProgram, fragShader);
        glLinkProgram(blitProgram);

        GLint success;
        glGetProgramiv(blitProgram, GL_LINK_STATUS, &success);
        if (!success) {
            char log[2048];
            glGetProgramInfoLog(blitProgram, 2048, nullptr, log);
            std::cerr << "Blit program link error:\n" << log << std::endl;
        }

        glDeleteShader(vertShader);
        glDeleteShader(fragShader);
        return success != 0;
    }

    void createFullscreenQuad() {
        float quadVertices[] = {
            // positions   // texCoords
            -1.0f,  1.0f,  0.0f, 1.0f,
            -1.0f, -1.0f,  0.0f, 0.0f,
             1.0f, -1.0f,  1.0f, 0.0f,

            -1.0f,  1.0f,  0.0f, 1.0f,
             1.0f, -1.0f,  1.0f, 0.0f,
             1.0f,  1.0f,  1.0f, 1.0f,
        };

        glGenVertexArrays(1, &quadVAO);
        glGenBuffers(1, &quadVBO);

        glBindVertexArray(quadVAO);
        glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);

        glBindVertexArray(0);
    }

    void createOutputTexture() {
        glGenTextures(1, &outputTexture);
        glBindTexture(GL_TEXTURE_2D, outputTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    }

    void createBuffers() {
        glGenBuffers(1, &simulationUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, simulationUBO);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(SimulationUBO), nullptr, GL_DYNAMIC_DRAW);

        glGenBuffers(1, &particleSSBO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, particleSSBO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, MAX_PARTICLES * sizeof(GPUParticle), nullptr, GL_DYNAMIC_DRAW);

        glGenBuffers(1, &flashSSBO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, flashSSBO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, MAX_FLASHES * sizeof(FusionFlash), nullptr, GL_DYNAMIC_DRAW);
    }
};

#endif // RAY_TRACING_H
