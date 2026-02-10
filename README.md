# TOKAMAK FUSION REACTOR SIMULATION (3D RAYTRACED)

A comprehensive **3D volumetric simulation** of a Tokamak nuclear fusion reactor, featuring realistic plasma physics, magnetic confinement, and a high-fidelity **ray-traced rendering engine** using OpenGL Compute Shaders.

![Project Status](https://img.shields.io/badge/status-active-brightgreen)
![OpenGL](https://img.shields.io/badge/OpenGL-4.6-blue)
![Language](https://img.shields.io/badge/C++-17-blue)

PROGRESS SO FAR:



https://github.com/user-attachments/assets/cf4fcc6d-587e-4783-8f35-de0766aaf0fe



## OVERVIEW

This project simulates the complex behavior of plasma within a tokamak reactor. Originally a 2D cross-sectional simulation, it has been upgraded to a full **3D volumetric experience**.

The simulation combines:
1.  **Plasma Physics Engine**: Calculates particle trajectories (Lorentz force), Coulomb collisions, and D-T fusion reactions.
2.  **Volumetric Raytracing**: A custom compute shader (`tokamak_raytrace.comp`) performs real-time raymarching to render the plasma density, temperature gradients, and the reactor vessel with physically-based shading.
3.  **Interactive Controls**: Full orbit camera control and real-time parameter tuning via ImGui.

## PROJECT STRUCTURE

```
project/
├── CMakeLists.txt              # Build configuration (CMake)
├── main.cpp                    # Application entry, windowing, and main loop
├── tokamak_raytrace.comp       # GLSL Compute Shader for volumetric raytracing
├── camera.h                    # Orbit camera implementation
├── plasma_physics.h            # core physics engine (particle dynamics, fusion)
├── magnetic_field.h            # Magnetic field calculations (toroidal + poloidal)
├── tokamak_geometry.h          # Reactor geometry definition
├── particle.h                  # Particle data structures
├── ray_tracing.cpp             # Host-side raytracing setup and buffer management
├── particle.vert / .frag       # (Legacy) Basic rasterization shaders
├── external/                   # Dependencies (GLFW, GLAD, GLM)
├── imgui-master/               # ImGui library for UI
└── PHYSICS_LECTURE.md          # Comprehensive physics documentation
```

## FEATURES

### Physics Simulation
- **3D Magnetic Confinement**: Particles follow helical paths along the two magnetic field lines (toroidal and poloidal (calculated)) in a torus.
- **Fusion Reactions**: Deuterium-Tritium (D-T) fusion producing Alpha particles and Neutrons.
- **Thermal Plasma**: Maxwellian velocity distribution at 150 million Kelvin.
- **Particle Behavior**: Coulomb repulsion, Debye shielding, and boundary reflection.

### Visual & Rendering
- **Volumetric Raytracing**: Real-time accumulation of plasma density and emission using raymarching.
- **Dynamic Lighting**: Plasma glows based on density and temperature; fusion events create bright flashes.
- **Reactor Vessel**: Physically-based rendering (PBR) of the tokamak interior.
- **ImGui Integration**: Real-time control over simulation parameters (under development).

###  Camera & Controls
- **Orbit Camera**:
    - **Rotate**: Left Mouse Button + Drag
    - **Pan**: Right Mouse Button + Drag
    - **Zoom**: Mouse Scroll Wheel
- **Keyboard**:
    - **ESC**: Exit simulation

## COMPILATION

The project now uses **CMake** for building.

### Prerequisites
- **OpenGL 4.3+** (Required for Compute Shaders)
- **CMake 3.20+**
- **C++17 Compiler** (MSVC recommended on Windows)

### Build Instructions (Windows)

1.  Open the project folder in a terminal.
2.  Create a build directory:
    ```bash
    mkdir build
    cd build
    ```
3.  Configure with CMake:
    ```bash
    cmake ..
    ```
4.  Build the project:
    ```bash
    cmake --build . --config Release
    ```
5.  Run the executable:
    ```bash
    ./Release/FusionTokamakSim.exe
    ```

*Note: The current `CMakeLists.txt` is configured with prebuilt libraries for Windows. Linux/macOS users may need to adjust the library paths or use system-installed packages for GLFW/GLM.*

## PHYSICS DOCUMENTATION

See **PHYSICS_LECTURE.md** for:
- Detailed derivations of magnetic confinement.
- Fusion cross-sections and reaction rates.
- Plasma stability and tokamak geometry.

## SCIENTIFIC ACCURACY & LIMITATIONS

This simulation implements:
- ✓ Lorentz force (v × B)
- ✓ D-T Fusion energy release (17.6 MeV)
- ✓ Volumetric plasma density visualization
- ✓ Toroidal confinement geometry

**Simplifications**:
- The raytracer visualizes a density field derived from particle positions/functions, which is an approximation of the discrete particle data for performance.
- Relativistic effects are ignored.
- Full Magnetohydrodynamics (MHD) fluid simulation is approximated by particle kinetics.

## LICENSE

This code is provided for educational and research purposes.

---

**Note**: Rendering volumetric plasma is GPU-intensive. If you experience low framerates, try reducing the window size or adjusting the raymarching step size in `tokamak_raytrace.comp`.
