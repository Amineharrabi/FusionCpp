#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>

struct OrbitCamera {
    float yaw   = 0.0f;      
    float pitch  = 0.4f;      
    float distance = 4.0f;    

    glm::vec3 target = glm::vec3(0.0f, 0.0f, 0.0f);

    glm::vec3 panOffset = glm::vec3(0.0f);

    float minDistance = 1.0f;
    float maxDistance = 15.0f;
    float minPitch = -1.5f;   
    float maxPitch = 1.5f;    

    float targetYaw = 0.0f;
    float targetPitch = 0.4f;
    float targetDistance = 4.0f;
    glm::vec3 targetPanOffset = glm::vec3(0.0f);
    float smoothFactor = 8.0f;

    bool leftDragging = false;
    bool rightDragging = false;
    double lastMouseX = 0.0;
    double lastMouseY = 0.0;
    float orbitSensitivity = 0.005f;
    float panSensitivity = 0.005f;
    float zoomSensitivity = 0.3f;

    float fov = 45.0f;
    float nearPlane = 0.01f;
    float farPlane = 100.0f;

    glm::vec3 getPosition() const {
        float x = distance * std::cos(pitch) * std::sin(yaw);
        float y = distance * std::sin(pitch);
        float z = distance * std::cos(pitch) * std::cos(yaw);
        return target + panOffset + glm::vec3(x, y, z);
    }

    glm::mat4 getViewMatrix() const {
        glm::vec3 pos = getPosition();
        glm::vec3 center = target + panOffset;
        return glm::lookAt(pos, center, glm::vec3(0.0f, 1.0f, 0.0f));
    }

    glm::mat4 getProjectionMatrix(float aspect) const {
        return glm::perspective(glm::radians(fov), aspect, nearPlane, farPlane);
    }

    glm::mat4 getInverseViewProjection(float aspect) const {
        glm::mat4 vp = getProjectionMatrix(aspect) * getViewMatrix();
        return glm::inverse(vp);
    }

    void update(float dt) {
        float t = 1.0f - std::exp(-smoothFactor * dt);
        yaw += (targetYaw - yaw) * t;
        pitch += (targetPitch - pitch) * t;
        distance += (targetDistance - distance) * t;
        panOffset += (targetPanOffset - panOffset) * t;
    }

    void onMouseButton(int button, int action, double mouseX, double mouseY) {
        if (button == 0) { // Left
            leftDragging = (action == 1);
            lastMouseX = mouseX;
            lastMouseY = mouseY;
        }
        if (button == 1) { // Right
            rightDragging = (action == 1);
            lastMouseX = mouseX;
            lastMouseY = mouseY;
        }
    }

    void onMouseMove(double mouseX, double mouseY) {
        float dx = (float)(mouseX - lastMouseX);
        float dy = (float)(mouseY - lastMouseY);
        lastMouseX = mouseX;
        lastMouseY = mouseY;

        if (leftDragging) {
            targetYaw -= dx * orbitSensitivity;
            targetPitch += dy * orbitSensitivity;
            if (targetPitch < minPitch) targetPitch = minPitch;
            if (targetPitch > maxPitch) targetPitch = maxPitch;
        }

        if (rightDragging) {
            glm::vec3 forward = glm::normalize(target + panOffset - getPosition());
            glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
            glm::vec3 up = glm::normalize(glm::cross(right, forward));
            targetPanOffset += (-right * dx + up * dy) * panSensitivity * distance;
        }
    }

    void onScroll(double yoffset) {
        targetDistance -= (float)yoffset * zoomSensitivity;
        if (targetDistance < minDistance) targetDistance = minDistance;
        if (targetDistance > maxDistance) targetDistance = maxDistance;
    }
};

#endif // CAMERA_H
