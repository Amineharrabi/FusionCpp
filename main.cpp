#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

// Include all simulation modules
#include "particle.h"
#include "tokamak_geometry.h"
#include "magnetic_field.h"
#include "plasma_physics.h"

// Shader loading utility
std::string loadShaderSource(const char *path)
{
    std::ifstream file(path);
    std::stringstream buf;
    buf << file.rdbuf();
    return buf.str();
}

GLuint compileShader(GLenum type, const char *src)
{
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cerr << "Shader compile error: " << infoLog << std::endl;
    }
    return shader;
}

int main()
{
    // ================== GLFW/OpenGL INITIALIZATION ==================
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow *window = glfwCreateWindow(1200, 800, "Tokamak Fusion Reactor Simulation", nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    int windowWidth = 1200, windowHeight = 800;
    glViewport(0, 0, windowWidth, windowHeight);

    // ================== SHADER COMPILATION ==================
    std::string vertSrc = loadShaderSource("particle.vert");
    std::string fragSrc = loadShaderSource("particle.frag");
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertSrc.c_str());
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragSrc.c_str());
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // ================== TOKAMAK SETUP ==================
    std::cout << "\n============================================" << std::endl;
    std::cout << "TOKAMAK FUSION REACTOR SIMULATION" << std::endl;
    std::cout << "============================================\n" << std::endl;
    
    // Initialize geometry (reactor structure)
    TokamakGeometry tokamak;
    std::cout << "Tokamak geometry initialized:" << std::endl;
    std::cout << "  Major radius: " << tokamak.majorRadius << " m (simulation units)" << std::endl;
    std::cout << "  Minor radius: " << tokamak.minorRadius << " m" << std::endl;
    std::cout << "  Plasma elongation: " << tokamak.plasmaElongation << std::endl;
    std::cout << "  Triangularity: " << tokamak.plasmaTriangularity << std::endl;
    
    // Initialize magnetic field
    MagneticField magneticField(tokamak.majorRadius, tokamak.minorRadius, 8.0f);
    std::cout << "\nMagnetic field configured:" << std::endl;
    std::cout << "  Toroidal field: " << magneticField.B_toroidal << " T" << std::endl;
    std::cout << "  Poloidal field: " << magneticField.B_poloidal << " T" << std::endl;
    std::cout << "  Safety factor q: " << magneticField.safetyFactor << std::endl;
    
    PlasmaPhysics plasmaPhysics(magneticField, tokamak);
    std::cout << "\nready" << std::endl;
    
    int numDeuterium = 400;  
    int numTritium = 400;
    
    std::vector<Particle> particles = plasmaPhysics.createThermalPlasma(numDeuterium, numTritium);
    
    std::cout << "\nInitial plasma created:" << std::endl;
    std::cout << "  Deuterium ions: " << numDeuterium << std::endl;
    std::cout << "  Tritium ions: " << numTritium << std::endl;
    std::cout << "  Total particles: " << particles.size() << std::endl;
    std::cout << "\nfusion events yellow flashes" << std::endl;
    
    // ================== RENDERING SETUP ==================
    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    // Timing
    double lastTime = glfwGetTime();
    int fusionCount = 0;
    double lastFusionTime = lastTime;

    bool simulationRunning = false;
    float injectionKick = 0.15f;
    std::mt19937 uiRng(std::random_device{}());

    // ================== MAIN SIMULATION LOOP ==================
    while (!glfwWindowShouldClose(window))
    {
        double currentTime = glfwGetTime();
        float deltaTime = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;
        
        // Limit dt for numerical stability
        if (deltaTime > 0.033f) deltaTime = 0.033f;

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Plasma Controls");
        float timeScale = plasmaPhysics.getTimeScale();
        float plasmaTemperature = plasmaPhysics.getPlasmaTemperature();
        float particleDensity = plasmaPhysics.getParticleDensity();
        float velocityScale = plasmaPhysics.getVelocityScale();
        float fusionBoost = plasmaPhysics.getFusionBoost();
        float maxFusionFractionPerStep = plasmaPhysics.getMaxFusionFractionPerStep();
        float confinementStrength = plasmaPhysics.getConfinementStrength();
        float coreAttractionStrength = plasmaPhysics.getCoreAttractionStrength();
        float driftOmega = plasmaPhysics.getDriftOmega();
        float wallLossProbability = plasmaPhysics.getWallLossProbability();
        bool enableCoulomb = plasmaPhysics.getEnableCoulomb();

        if (!simulationRunning) {
            ImGui::Text("Status: Paused");
            ImGui::SliderFloat("Injection Kick", &injectionKick, 0.0f, 2.0f, "%.3f");
            if (ImGui::Button("Start Injection")) {
                std::uniform_real_distribution<float> angleDist(0.0f, 2.0f * 3.1415926535f);
                for (auto& p : particles) {
                    if (!p.active) continue;
                    if (p.type != Particle::DEUTERIUM && p.type != Particle::TRITIUM) continue;
                    float a = angleDist(uiRng);
                    p.vx += injectionKick * std::cos(a);
                    p.vy += injectionKick * std::sin(a);
                }
                simulationRunning = true;
            }
        } else {
            ImGui::Text("Status: Running");
        }

        if (ImGui::SliderFloat("Time Scale", &timeScale, 1e-4f, 1.0f, "%.6f", ImGuiSliderFlags_Logarithmic)) {
            plasmaPhysics.setTimeScale(timeScale);
        }
        if (ImGui::SliderFloat("Temperature (K)", &plasmaTemperature, 1e7f, 5e9f, "%.3e", ImGuiSliderFlags_Logarithmic)) {
            plasmaPhysics.setPlasmaTemperature(plasmaTemperature);
        }
        if (ImGui::SliderFloat("Density (m^-3)", &particleDensity, 1e18f, 1e22f, "%.3e", ImGuiSliderFlags_Logarithmic)) {
            plasmaPhysics.setParticleDensity(particleDensity);
        }
        if (ImGui::SliderFloat("Velocity Scale", &velocityScale, 1e-9f, 1e-4f, "%.3e", ImGuiSliderFlags_Logarithmic)) {
            plasmaPhysics.setVelocityScale(velocityScale);
        }
        if (ImGui::SliderFloat("Fusion Boost", &fusionBoost, 1.0f, 1e9f, "%.3e", ImGuiSliderFlags_Logarithmic)) {
            plasmaPhysics.setFusionBoost(fusionBoost);
        }
        if (ImGui::SliderFloat("Max Fusion Fraction/Step", &maxFusionFractionPerStep, 0.0f, 0.2f, "%.3f")) {
            plasmaPhysics.setMaxFusionFractionPerStep(maxFusionFractionPerStep);
        }

        if (ImGui::SliderFloat("Confinement Strength", &confinementStrength, 0.0f, 500.0f, "%.3f")) {
            plasmaPhysics.setConfinementStrength(confinementStrength);
        }

        if (ImGui::SliderFloat("Core Attraction", &coreAttractionStrength, 0.0f, 50.0f, "%.3f")) {
            plasmaPhysics.setCoreAttractionStrength(coreAttractionStrength);
        }
        if (ImGui::SliderFloat("Drift Omega", &driftOmega, 0.0f, 20.0f, "%.3f")) {
            plasmaPhysics.setDriftOmega(driftOmega);
        }
        if (ImGui::SliderFloat("Wall Loss Probability", &wallLossProbability, 0.0f, 1.0f, "%.3f")) {
            plasmaPhysics.setWallLossProbability(wallLossProbability);
        }
        if (ImGui::Checkbox("Enable Coulomb (slow)", &enableCoulomb)) {
            plasmaPhysics.setEnableCoulomb(enableCoulomb);
        }

        ImGui::Text("Particles: %d", (int)particles.size());
        ImGui::Text("Fusions: %d", fusionCount);
        ImGui::End();

        // ============ PHYSICS UPDATE ============
        if (simulationRunning) {
            plasmaPhysics.updateParticles(particles, deltaTime);
        }
        
        // Count active particles and fusion products
        int activeD = 0, activeT = 0, heliumCount = 0, neutronCount = 0;
        for (const auto& p : particles) {
            if (!p.active) continue;
            if (p.type == Particle::DEUTERIUM) activeD++;
            else if (p.type == Particle::TRITIUM) activeT++;
            else if (p.type == Particle::HELIUM) heliumCount++;
            else if (p.type == Particle::NEUTRON) neutronCount++;
        }
        
        int currentFusions = heliumCount;
        if (currentFusions > fusionCount) {
            std::cout << "FUSION occured Total fusions: " << currentFusions 
                     << " | D: " << activeD << " T: " << activeT 
                     << " | He: " << heliumCount << " n: " << neutronCount << std::endl;
            fusionCount = currentFusions;
            lastFusionTime = currentTime;
        }

        // ============ RENDERING ============
        glClearColor(0.02f, 0.02f, 0.05f, 1.0f); // Dark background
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        
        // Draw Tokamak structure
        tokamak.render(shaderProgram);
        
        // Draw particles
        glBindVertexArray(VAO);
        GLint colorLoc = glGetUniformLocation(shaderProgram, "particleColor");
        
        for (const auto& p : particles) {
            if (!p.active) continue;
            
            // Generate circle vertices for particle
            std::vector<float> verts = generateCircleVertices(p.x, p.y, p.radius, 16);
            
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), 
                        verts.data(), GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
            glEnableVertexAttribArray(0);
            
            // Add glow effect for recent fusion products
            float glowFactor = 1.0f;
            if (p.type == Particle::HELIUM && (currentTime - lastFusionTime) < 0.5f) {
                glowFactor = 2.0f; // Bright flash for new helium
            }
            
            glUniform4f(colorLoc, p.r * glowFactor, p.g * glowFactor, p.b * glowFactor, p.a);
            glDrawArrays(GL_TRIANGLE_FAN, 0, verts.size() / 2);
        }

        glBindVertexArray(0);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        // Exit on ESC key
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, true);
        }
    }

    std::cout << "\nSimulation ended." << std::endl;
    std::cout << "Final statistics:" << std::endl;
    std::cout << "  Total fusion reactions: " << fusionCount << std::endl;
    std::cout << "  Final particle count: " << particles.size() << std::endl;

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);
    tokamak.cleanup();
    glfwDestroyWindow(window);
    glfwTerminate();
    
    return 0;
}
