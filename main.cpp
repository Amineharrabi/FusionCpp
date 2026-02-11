#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "particle.h"
#include "tokamak_geometry.h"
#include "magnetic_field.h"
#include "plasma_physics.h"
#include "camera.h"
#include "ray_tracing.cpp"  
OrbitCamera g_camera;
int g_windowWidth = 1200;
int g_windowHeight = 800;

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) return;
    
    double mx, my;
    glfwGetCursorPos(window, &mx, &my);
    g_camera.onMouseButton(button, action, mx, my);
}

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) return;
    
    g_camera.onMouseMove(xpos, ypos);
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) return;
    
    g_camera.onScroll(yoffset);
}

void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
    g_windowWidth = width;
    g_windowHeight = height;
    glViewport(0, 0, width, height);
}

void fatalError(const char* msg) {
    std::cerr << "FATAL ERROR: " << msg << std::endl;
    std::cout << "Press Enter to exit..." << std::endl;
    std::cin.get();
    exit(-1);
}

int main()
{
    if (!glfwInit()) {
        fatalError("Failed to initialize GLFW");
    }
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(g_windowWidth, g_windowHeight, 
                                           "Tokamak Fusion Reactor — 3D Ray Tracing", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        fatalError("Failed to create GLFW window (OpenGL 4.3 required)");
    }
    glfwMakeContextCurrent(window);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        fatalError("Failed to initialize GLAD");
    }

    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GPU: " << glGetString(GL_RENDERER) << std::endl;

    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glViewport(0, 0, g_windowWidth, g_windowHeight);

    GPURayTracer rayTracer;
    if (!rayTracer.initialize(g_windowWidth, g_windowHeight)) {
        fatalError("Failed to initialize GPU ray tracer (check console for shader errors)");
    }

    std::cout << "\n============================================" << std::endl;
    std::cout << "TOKAMAK FUSION REACTOR — 3D SIMULATION" << std::endl;
    std::cout << "============================================\n" << std::endl;

    TokamakGeometry tokamak;
    std::cout << "Torus geometry:" << std::endl;
    std::cout << "  Major radius: " << tokamak.torusMajorR << std::endl;
    std::cout << "  Minor radius: " << tokamak.torusMinorR << std::endl;

    MagneticField magneticField(tokamak.torusMajorR, tokamak.torusMinorR, 8.0f);
    std::cout << "Magnetic field: Bt=" << magneticField.B_toroidal 
              << " T, Bp=" << magneticField.B_poloidal << " T" << std::endl;

    PlasmaPhysics plasmaPhysics(magneticField, tokamak);
    
    int numDeuterium = 4200;
    int numTritium = 4200;
    std::vector<Particle> particles = plasmaPhysics.createThermalPlasma(numDeuterium, numTritium);
    
    std::cout << "Initial plasma: " << numDeuterium << " D + " << numTritium << " T = " 
              << particles.size() << " particles" << std::endl;
    std::cout << "\nControls: LMB drag = orbit, Scroll = zoom, RMB drag = pan" << std::endl;
    std::cout << "Press Start Injection to begin fusion!" << std::endl;

    double lastTime = glfwGetTime();
    int fusionCount = 0;
    double lastFusionTime = lastTime;

    bool simulationRunning = false;
    float injectionKick = 0.25f;
    std::mt19937 uiRng(std::random_device{}());

    std::vector<FusionFlash> activeFlashes;
    const float flashDuration = 2.5f; 

    bool autoFuel = true;
    int fuelThreshold = 5000; 
    int fuelBatchSize = 1000;  
    float fuelCooldown = 0.0f;
    float fuelCooldownTime = 0.6f; 
    while (!glfwWindowShouldClose(window))
    {
        double currentTime = glfwGetTime();
        float deltaTime = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;
        if (deltaTime > 0.033f) deltaTime = 0.033f;

        g_camera.update(deltaTime);

        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        if (fbWidth != rayTracer.width || fbHeight != rayTracer.height) {
            if (fbWidth > 0 && fbHeight > 0) {
                rayTracer.resize(fbWidth, fbHeight);
                g_windowWidth = fbWidth;
                g_windowHeight = fbHeight;
            }
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Plasma Controls");

        float timeScale = plasmaPhysics.getTimeScale();
        float plasmaTemperature = plasmaPhysics.getPlasmaTemperature();
        float particleDensity = plasmaPhysics.getParticleDensity();
        float velocityScale = plasmaPhysics.getVelocityScale();
        float fusionBoost = plasmaPhysics.getFusionBoost();
        float maxFusionFrac = plasmaPhysics.getMaxFusionFractionPerStep();
        float confinement = plasmaPhysics.getConfinementStrength();
        float coreAttraction = plasmaPhysics.getCoreAttractionStrength();
        float driftOmega = plasmaPhysics.getDriftOmega();
        float wallLoss = plasmaPhysics.getWallLossProbability();
        bool coulomb = plasmaPhysics.getEnableCoulomb();

        if (!simulationRunning) {
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "Status: PAUSED");
            ImGui::SliderFloat("Injection Kick", &injectionKick, 0.0f, 2.0f, "%.3f");
            if (ImGui::Button("Start Injection")) {
                std::uniform_real_distribution<float> angleDist(0.0f, 2.0f * 3.14159f);
                std::uniform_real_distribution<float> angleDist2(0.0f, 2.0f * 3.14159f);
                for (auto& p : particles) {
                    if (!p.active) continue;
                    if (p.type != Particle::DEUTERIUM && p.type != Particle::TRITIUM) continue;
                    float a = angleDist(uiRng);
                    float b = angleDist2(uiRng);
                    float R = std::sqrt(p.x * p.x + p.z * p.z);
                    if (R > 1e-6f) {
                        p.vx += injectionKick * (-p.z / R);
                        p.vz += injectionKick * (p.x / R);
                    }
                    p.vy += injectionKick * 0.3f * std::sin(b);
                }
                simulationRunning = true;
            }
        } else {
            ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.5f, 1.0f), "Status: RUNNING");
        }

        ImGui::Separator();
        ImGui::Text("--- Physics ---");
        
        if (ImGui::SliderFloat("Time Scale", &timeScale, 1e-4f, 1.0f, "%.6f", ImGuiSliderFlags_Logarithmic))
            plasmaPhysics.setTimeScale(timeScale);
        if (ImGui::SliderFloat("Temperature (K)", &plasmaTemperature, 1e7f, 5e9f, "%.3e", ImGuiSliderFlags_Logarithmic))
            plasmaPhysics.setPlasmaTemperature(plasmaTemperature);
        if (ImGui::SliderFloat("Fusion Boost", &fusionBoost, 1.0f, 1e9f, "%.3e", ImGuiSliderFlags_Logarithmic))
            plasmaPhysics.setFusionBoost(fusionBoost);
        if (ImGui::SliderFloat("Confinement", &confinement, 0.0f, 500.0f, "%.1f"))
            plasmaPhysics.setConfinementStrength(confinement);
        if (ImGui::SliderFloat("Core Attraction", &coreAttraction, 0.0f, 50.0f, "%.1f"))
            plasmaPhysics.setCoreAttractionStrength(coreAttraction);
        if (ImGui::SliderFloat("Drift Omega", &driftOmega, 0.0f, 20.0f, "%.1f"))
            plasmaPhysics.setDriftOmega(driftOmega);

        ImGui::Separator();
        ImGui::Text("--- Torus Rendering ---");
        ImGui::SliderFloat("Torus Opacity", &tokamak.torusOpacity, 0.0f, 1.0f, "%.2f");

        ImGui::Separator();
        ImGui::Text("--- Fueling ---");
        ImGui::Checkbox("Auto-Fuel", &autoFuel);
        ImGui::SliderInt("Fuel Threshold", &fuelThreshold, 10, 5000);
        ImGui::SliderInt("Fuel Batch Size", &fuelBatchSize, 10, 1000);
        if (ImGui::Button("Manual Refuel")) {
            plasmaPhysics.injectFuel(particles, fuelBatchSize, fuelBatchSize);
            std::cout << "REFUELED: +" << fuelBatchSize << " D + " << fuelBatchSize << " T" << std::endl;
        }

        int activeD = 0, activeT = 0, heliumCount = 0, neutronCount = 0, totalActive = 0;
        for (const auto& p : particles) {
            if (!p.active) continue;
            totalActive++;
            if (p.type == Particle::DEUTERIUM) activeD++;
            else if (p.type == Particle::TRITIUM) activeT++;
            else if (p.type == Particle::HELIUM) heliumCount++;
            else if (p.type == Particle::NEUTRON) neutronCount++;
        }

        ImGui::Separator();
        ImGui::Text("--- Statistics ---");
        ImGui::Text("Active particles: %d", totalActive);
        ImGui::TextColored(ImVec4(0.3f, 0.6f, 1.0f, 1.0f), "  Deuterium: %d", activeD);
        ImGui::TextColored(ImVec4(0.6f, 0.3f, 1.0f, 1.0f), "  Tritium: %d", activeT);
        ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.3f, 1.0f), "  Helium-4: %d", heliumCount);
        ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "  Neutrons: %d", neutronCount);
        ImGui::Text("Fusion events: %d", fusionCount);
        ImGui::Text("Active flashes: %d", (int)activeFlashes.size());
        ImGui::Text("FPS: %.1f", 1.0f / deltaTime);

        ImGui::End();

        if (simulationRunning) {
            int preFusionHe = heliumCount;

            plasmaPhysics.updateParticles(particles, deltaTime);

            int postFusionHe = 0;
            for (const auto& p : particles) {
                if (p.active && p.type == Particle::HELIUM) postFusionHe++;
            }

            int newFusions = postFusionHe - preFusionHe;
            if (newFusions > 0) {
                fusionCount += newFusions;
                lastFusionTime = currentTime;

                for (int i = (int)particles.size() - 1; i >= 0 && newFusions > 0; --i) {
                    if (particles[i].active && particles[i].type == Particle::HELIUM) {
                        FusionFlash flash;
                        flash.px = particles[i].x;
                        flash.py = particles[i].y;
                        flash.pz = particles[i].z;
                        flash.age = 0.0f;
                        flash.r = 1.0f;
                        flash.g = 0.95f;
                        flash.b = 0.4f;
                        flash.intensity = 2.0f;
                        activeFlashes.push_back(flash);
                        newFusions--;
                    }
                }

                std::cout << "fusion happned Total: " << fusionCount
                          << " Deuterium:" << activeD << " T:" << activeT
                          << " Helium:" << postFusionHe << std::endl;
            }

            if (autoFuel) {
                fuelCooldown -= deltaTime;
                int curD = 0, curT = 0;
                for (const auto& p : particles) {
                    if (!p.active) continue;
                    if (p.type == Particle::DEUTERIUM) curD++;
                    else if (p.type == Particle::TRITIUM) curT++;
                }
                if ((curD < fuelThreshold || curT < fuelThreshold) && fuelCooldown <= 0.0f) {
                    plasmaPhysics.injectFuel(particles, fuelBatchSize, fuelBatchSize);
                    fuelCooldown = fuelCooldownTime;
                    std::cout << "autoFuel: +" << fuelBatchSize << " D + " << fuelBatchSize 
                              << " T (D was " << curD << ", T was " << curT << ")" << std::endl;
                }
            }

            if (particles.size() > 15000) {
                particles.erase(
                    std::remove_if(particles.begin(), particles.end(),
                                   [](const Particle& p) { return !p.active; }),
                    particles.end()
                );
            }
        }

        for (auto& flash : activeFlashes) {
            flash.age += deltaTime / flashDuration;
        }
        activeFlashes.erase(
            std::remove_if(activeFlashes.begin(), activeFlashes.end(),
                           [](const FusionFlash& f) { return f.age >= 1.0f; }),
            activeFlashes.end()
        );


        std::vector<GPUParticle> gpuParticles; 


        std::vector<FusionFlash> gpuFlashes = activeFlashes;
        if ((int)gpuFlashes.size() > GPURayTracer::MAX_FLASHES) {
            gpuFlashes.resize(GPURayTracer::MAX_FLASHES);
        }

        float aspect = (float)g_windowWidth / (float)g_windowHeight;
        glm::mat4 invVP = g_camera.getInverseViewProjection(aspect);
        glm::vec3 camPos = g_camera.getPosition();

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        rayTracer.render(
            invVP,
            camPos,
            tokamak.torusMajorR,
            tokamak.torusMinorR,
            tokamak.torusOpacity,
            (float)currentTime,
            gpuParticles,
            gpuFlashes,
            totalActive
        );

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, true);
        }
    }

    std::cout << "\nSimulation ended." << std::endl;
    std::cout << "Total fusion reactions: " << fusionCount << std::endl;
    std::cout << "Final particle count: " << particles.size() << std::endl;

    rayTracer.cleanup();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
