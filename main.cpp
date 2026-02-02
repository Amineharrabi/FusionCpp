
#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include "particle.h"

int main()
{
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
    GLFWwindow *window = glfwCreateWindow(800, 600, "Fusion Reactor Simulation", nullptr, nullptr);
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

    // Define a particle and generate its circle vertices
    Particle particle = {0.0f, 0.0f, 0.1f, 1.0f, 0.5f, 0.2f, 1.0f}; // Center in NDC space
    std::vector<float> circleVertices = generateCircleVertices(particle.x, particle.y, particle.radius, 64);

    // --- Shader loading utility ---
    auto loadShaderSource = [](const char *path) -> std::string
    {
        std::ifstream file(path);
        std::stringstream buf;
        buf << file.rdbuf();
        return buf.str();
    };

    auto compileShader = [](GLenum type, const char *src) -> GLuint
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
    };

    // Load and compile shaders
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

    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, circleVertices.size() * sizeof(float), circleVertices.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Animation variables (optional, for moving the particle)
    float speed = 2.0f; // Units per second in NDC space
    double lastTime = glfwGetTime();

    glPointSize(10.0f);
    while (!glfwWindowShouldClose(window))
    {
        double currentTime = glfwGetTime();
        float deltaTime = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;

        // Update particle position (move horizontally)
        particle.x += speed * deltaTime;

        // Wrap around screen
        if (particle.x > 1.0f)
            particle.x = -1.0f;

        // Regenerate vertices with new position
        circleVertices = generateCircleVertices(particle.x, particle.y, particle.radius, 64);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, circleVertices.size() * sizeof(float), circleVertices.data());
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);
        glDrawArrays(GL_POINTS, 0, 1);
        glBindVertexArray(0);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
