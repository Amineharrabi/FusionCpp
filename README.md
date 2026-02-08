# TOKAMAK FUSION REACTOR SIMULATION

A comprehensive 2D cross-sectional simulation of a Tokamak nuclear fusion reactor with realistic plasma physics and magnetic confinement.

## PROJECT STRUCTURE

```
project/
├── main.cpp                    # Main simulation loop and rendering
├── particle.h                  # Particle structure and physics constants
├── tokamak_geometry.h          # Reactor geometry and D-shape plasma
├── magnetic_field.h            # Magnetic field calculations (toroidal + poloidal)
├── plasma_physics.h            # Particle dynamics, collisions, and fusion
├── particle.vert              # Vertex shader (simple passthrough)
├── particle.frag              # Fragment shader (simple coloring)
├── PHYSICS_LECTURE.md         # Comprehensive physics documentation
└── README.md                  # This file
```

## FEATURES

### Physics Simulation

- **Magnetic Confinement**: Realistic Lorentz force from combined toroidal and poloidal fields
- **Coulomb Interactions**: Particle-particle repulsion with Debye shielding
- **Fusion Reactions**: D-T fusion producing alpha particles and neutrons with proper energy distribution
- **Thermal Plasma**: Maxwellian velocity distribution at 150 million Kelvin
- **Boundary Conditions**: Particle reflection at plasma edge with energy loss

### Visual Features

- **Tokamak Geometry**: D-shaped plasma cross-section with elongation and triangularity
- **Magnetic Coils**: Toroidal field coils visualized around vessel
- **Divertor Region**: Lower exhaust region for particle removal
- **Real-time Particle Visualization**: Color-coded particles (blue=Deuterium, purple=Tritium, yellow=Helium, gray=Neutron)
- **Fusion Flash**: Visual burst when fusion reactions occur

## COMPILATION

### Prerequisites

- C++11 or later compiler (g++, clang++)
- GLFW 3.x
- GLAD (OpenGL loader)
- OpenGL 3.3+

### Linux/MacOS

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt-get install libglfw3-dev

# Compile
g++ -std=c++11 main.cpp -o tokamak_sim \
    -lglfw -lGL -ldl -lpthread \
    -I. -O2

# Run
./tokamak_sim
```

### MacOS Specific

```bash
# Install GLFW via Homebrew
brew install glfw

# Compile
g++ -std=c++11 main.cpp -o tokamak_sim \
    -I/opt/homebrew/include -L/opt/homebrew/lib \
    -lglfw -framework OpenGL -framework Cocoa \
    -framework IOKit -framework CoreVideo \
    -I. -O2

# Run
./tokamak_sim
```

### Windows (MinGW)

```bash
g++ -std=c++11 main.cpp -o tokamak_sim.exe -lglfw3 -lopengl32 -lgdi32 -I. -O2
```

## SHADER FILES

You need two simple shader files in the same directory:

**particle.vert**:

```glsl
#version 330 core
layout (location = 0) in vec2 aPos;
void main()
{
    gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);
}
```

**particle.frag**:

```glsl
#version 330 core
out vec4 FragColor;
uniform vec4 particleColor;
void main()
{
    FragColor = particleColor;
}
```

## USAGE

### Running the Simulation

```bash
./tokamak_sim
```

### Controls

- **ESC**: Exit simulation

### Console Output

The simulation prints:

- Tokamak parameters (radii, field strengths)
- Initial particle counts
- Fusion events in real-time with particle statistics

Example output:

```
============================================
TOKAMAK FUSION REACTOR SIMULATION
============================================

Tokamak geometry initialized:
  Major radius: 0.6 m (simulation units)
  Minor radius: 0.25 m
  Plasma elongation: 1.7
  Triangularity: 0.33

Magnetic field configured:
  Toroidal field: 8 T
  Poloidal field: 0.111111 T
  Safety factor q: 3

Plasma physics engine ready.

Initial plasma created:
  Deuterium ions: 100
  Tritium ions: 100
  Total particles: 200

Simulation started. Watch for fusion events (yellow flashes)!
============================================

FUSION EVENT! Total fusions: 1 | D: 98 T: 99 | He: 1 n: 1
FUSION EVENT! Total fusions: 2 | D: 96 T: 98 | He: 2 n: 2
...
```

## TUNING PARAMETERS

### In main.cpp

**Particle Count** (Line 110):

```cpp
int numDeuterium = 100;  // Increase for more particles (100-500 recommended)
int numTritium = 100;
```

**Time Scale** (plasma_physics.h, Line 29):

```cpp
timeScale(1e-9f)  // Increase for faster dynamics, decrease for stability
```

### In tokamak_geometry.h

**Reactor Size**:

```cpp
majorRadius(0.6f)      // Increase for larger reactor
minorRadius(0.25f)
plasmaElongation(1.7f)  // 1.0 = circular, >1.0 = elongated
plasmaTriangularity(0.33f)  // D-shape parameter
```

### In magnetic_field.h

**Field Strength**:

```cpp
MagneticField magneticField(R, a, 8.0f);  // Last parameter = toroidal field (Tesla)
```

**Safety Factor** (Line 33):

```cpp
safetyFactor(3.0f)  // Typical: 2-4, higher = more stable
```

### In plasma_physics.h

**Plasma Temperature** (Line 30):

```cpp
plasmaTemperature(1.5e8f)  // Kelvin (1.5e8 = 150 million K ≈ 13 keV)
```

**Fusion Probability Scaling** (Line 261):

```cpp
float fusionChance = crossSection * 1e28f;  // Increase multiplier for more fusions
```

## PHYSICS DOCUMENTATION

See **PHYSICS_LECTURE.md** for:

- Complete derivations of all equations
- Physical constants and typical values
- Detailed explanation of magnetic confinement
- Fusion cross-sections and reaction rates
- Numerical methods and stability

## VISUALIZATION GUIDE

### Color Coding

- **Cyan glow**: Plasma boundary
- **Dark gray**: Vacuum vessel
- **Orange**: Magnetic coils
- **Red**: Divertor region
- **Light blue particles**: Deuterium ions
- **Purple particles**: Tritium ions
- **Bright yellow particles**: Helium-4 (alpha particles from fusion)
- **Gray particles**: Neutrons (move unaffected by magnetic field)

### What to Look For

1. **Confinement**: Particles should gyrate in small circles and stay within the plasma boundary
2. **Boundary Collisions**: Particles hitting the edge should reflect back
3. **Fusion Events**: Watch for yellow flashes when D+T fuse into He+n
4. **Neutron Escape**: Neutrons (gray) will drift out since they're uncharged

## PERFORMANCE NOTES

- **Computational Complexity**: O(N²) for particle interactions
- **Recommended particle count**: 100-300 for smooth real-time
- **For 500+ particles**: Consider GPU acceleration or spatial hashing
- **Frame rate**: Should maintain 30+ FPS on modern hardware

## EXTENDING THE SIMULATION

### Adding Electrons

Uncomment electron creation in `plasma_physics.h`:

```cpp
// Add quasi-neutrality
for (int i = 0; i < numDeuterium + numTritium; ++i) {
    particles.push_back(createParticle(Particle::ELECTRON, x, y, vx, vy));
}
```

### Upgrading to 3D

1. Extend geometry to full torus (add toroidal angle)
2. Add third velocity component (vz)
3. Implement 3D rendering with proper camera
4. Consider instanced rendering for performance

### Adding Heating Systems

Implement:

- **Neutral Beam Injection**: High-energy neutral atoms
- **RF Heating**: Cyclotron resonance heating
- **Ohmic Heating**: From plasma current

### Advanced Physics

- **MHD Instabilities**: Kink modes, tearing modes
- **Turbulent Transport**: Anomalous diffusion
- **Alpha Heating**: Fusion products heating plasma
- **Bootstrap Current**: Self-generated current

## TROUBLESHOOTING

### Compilation Errors

- **GLFW not found**: Install GLFW development libraries
- **GLAD errors**: Make sure glad.h is in include path
- **Undefined references**: Check linker flags (-lglfw -lGL)

### Runtime Issues

- **Particles escape immediately**: Decrease timeScale in plasma_physics.h
- **No fusion events**: Increase particle count or fusion probability scaling
- **Crash/segfault**: Check OpenGL context creation and shader compilation
- **Black screen**: Verify shader files exist and compile correctly

### Physics Issues

- **Unrealistic behavior**: Check that physical constants match reality
- **Numerical instability**: Reduce time step or increase damping
- **Particles stuck**: Verify boundary reflection logic

## SCIENTIFIC ACCURACY

This simulation implements:
✓ Lorentz force (v × B)
✓ Coulomb repulsion with Debye shielding
✓ Magnetic mirror force
✓ D-T fusion with correct energy release
✓ Momentum and energy conservation in fusion
✓ Realistic Tokamak geometry
✓ Thermal velocity distributions

Simplifications made:

- 2D cross-section (not full 3D torus)
- Classical particle dynamics (no quantum effects)
- Simplified fusion cross-section
- No collective plasma effects (waves, instabilities)
- No radiation losses

## REFERENCES

1. Wesson, J. "Tokamaks" (4th Edition) - Standard tokamak textbook
2. Freidberg, J.P. "Plasma Physics and Fusion Energy" - Comprehensive fusion physics
3. Chen, F.F. "Introduction to Plasma Physics" - Plasma fundamentals
4. Stacey, W.M. "Fusion Plasma Physics" - Advanced plasma theory

## LICENSE

This code is provided for educational and research purposes.

## AUTHOR

Created as a comprehensive physics simulation demonstrating magnetic confinement fusion.

---

**Note**: This is a simplified educational simulation. Real fusion reactors involve vastly more complex physics including turbulence, MHD instabilities, radiation transport, and require supercomputer-scale simulations.

/\*\*

- PLASMA PHYSICS ENGINE
-
- This module handles:
- 1.  Particle dynamics under magnetic confinement
- 2.  Coulomb interactions between charged particles
- 3.  Fusion reactions when conditions are met
- 4.  Plasma boundary interactions
- 5.  Energy and momentum conservation
      \*/



