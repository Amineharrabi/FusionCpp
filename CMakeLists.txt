cmake_minimum_required(VERSION 3.20)

project(FusionTokamakSim LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(FusionTokamakSim
    main.cpp
    particle.cpp
    external/glad/src/glad.c

    imgui-master/imgui.cpp
    imgui-master/imgui_draw.cpp
    imgui-master/imgui_tables.cpp
    imgui-master/imgui_widgets.cpp
    imgui-master/imgui_demo.cpp
    imgui-master/backends/imgui_impl_glfw.cpp
    imgui-master/backends/imgui_impl_opengl3.cpp
)

target_include_directories(FusionTokamakSim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/glfw/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/glm/g-truc-glm-a532f5b
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui-master
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui-master/backends
)

target_compile_definitions(FusionTokamakSim PRIVATE IMGUI_IMPL_OPENGL_LOADER_GLAD)

# --- GLFW (prebuilt) ---
set(GLFW_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/glfw/lib")

if(WIN32)
    # Prefer dynamic import library (pairs with glfw3.dll)
    if(EXISTS "${GLFW_LIB_DIR}/glfw3dll.lib")
        set(GLFW_LIB_FILE "${GLFW_LIB_DIR}/glfw3dll.lib")
    elseif(EXISTS "${GLFW_LIB_DIR}/glfw3.lib")
        set(GLFW_LIB_FILE "${GLFW_LIB_DIR}/glfw3.lib")
    elseif(EXISTS "${GLFW_LIB_DIR}/glfw3_mt.lib")
        set(GLFW_LIB_FILE "${GLFW_LIB_DIR}/glfw3_mt.lib")
    else()
        message(FATAL_ERROR "Could not find a GLFW .lib in: ${GLFW_LIB_DIR}")
    endif()

    add_library(glfw_prebuilt STATIC IMPORTED GLOBAL)
    set_target_properties(glfw_prebuilt PROPERTIES
        IMPORTED_LOCATION "${GLFW_LIB_FILE}"
    )

    target_link_libraries(FusionTokamakSim PRIVATE
        glfw_prebuilt
        opengl32
    )

    # Copy runtime assets next to the .exe
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

    if(EXISTS "${GLFW_LIB_DIR}/glfw3.dll")
        add_custom_command(TARGET FusionTokamakSim POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${GLFW_LIB_DIR}/glfw3.dll"
                "$<TARGET_FILE_DIR:FusionTokamakSim>/glfw3.dll"
        )
    endif()
else()
    message(FATAL_ERROR "This CMakeLists currently supports the prebuilt Windows GLFW in external/glfw/lib.")
endif()
