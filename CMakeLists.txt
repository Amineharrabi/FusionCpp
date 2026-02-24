cmake_minimum_required(VERSION 3.20)

project(FusionTokamakSim LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
)

FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/Dav1dde/glad.git
    GIT_TAG v0.1.36
)

FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 0.9.9.8
)

FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.90.8
)

set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(glfw glad glm imgui)

add_executable(FusionTokamakSim
    main.cpp
    particle.cpp

    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)

target_include_directories(FusionTokamakSim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)

target_compile_definitions(FusionTokamakSim PRIVATE IMGUI_IMPL_OPENGL_LOADER_GLAD)

find_package(OpenGL REQUIRED)

target_link_libraries(FusionTokamakSim PRIVATE
    glfw
    glad
    OpenGL::GL
)

if(TARGET glm::glm)
    target_link_libraries(FusionTokamakSim PRIVATE glm::glm)
elseif(TARGET glm)
    target_link_libraries(FusionTokamakSim PRIVATE glm)
else()
    target_include_directories(FusionTokamakSim PRIVATE ${glm_SOURCE_DIR})
endif()

add_custom_command(TARGET FusionTokamakSim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/particle.vert"
        "$<TARGET_FILE_DIR:FusionTokamakSim>/particle.vert"
)
add_custom_command(TARGET FusionTokamakSim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/particle.frag"
        "$<TARGET_FILE_DIR:FusionTokamakSim>/particle.frag"
)
add_custom_command(TARGET FusionTokamakSim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/tokamak_raytrace.comp"
        "$<TARGET_FILE_DIR:FusionTokamakSim>/tokamak_raytrace.comp"
)
