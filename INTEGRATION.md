# INTEGRATION GUIDE: HOW TO ADD PARTICLES TO YOUR SIMULATION

## OVERVIEW

Your Tokamak simulation is now fully modular with separate files for geometry, physics, and magnetic fields. Here's exactly what each file does and where particles are created.

## FILE RESPONSIBILITIES

### 1. **particle.h** - Particle Definition
**What it does:**
- Defines the `Particle` structure with position, velocity, mass, charge
- Contains physical constants (masses, charges, Boltzmann constant, etc.)
- Provides `createParticle()` function to make particles of different types
- Includes `generateCircleVertices()` for rendering

**You don't need to modify this** - it's a library of particle utilities.

### 2. **tokamak_geometry.h** - Reactor Shape
**What it does:**
- Creates the D-shaped plasma boundary (cross-section)
- Draws vacuum vessel walls
- Positions magnetic coils
- Draws divertor region
- Provides `isInsidePlasma()` to check if a point is in the plasma

**Purpose:** Pure geometry and rendering - no particle physics here.

### 3. **magnetic_field.h** - Magnetic Field Calculations
**What it does:**
- Calculates toroidal field (B_φ) - the main confining field
- Calculates poloidal field (B_θ) - from plasma current
- Provides `getTotalField(x, y)` - returns B field vector at any point
- Calculates Larmor radius and magnetic pressure
- Implements `calculateLorentzForce()` - the F = q(v × B) force
- Implements `calculateMirrorForce()` - reflection from strong field regions

**Purpose:** All magnetic field physics - used by plasma_physics.h

### 4. **plasma_physics.h** - The Physics Engine (★ THIS IS WHERE PARTICLES ARE CREATED AND UPDATED ★)
**What it does:**
- **`updateParticles()`** - Main update loop that:
  - Applies magnetic forces (Lorentz + mirror)
  - Applies Coulomb forces between particles
  - Checks for and executes fusion reactions
  - Updates positions and velocities
  - Handles boundary collisions

- **`createThermalPlasma()`** - ★ CREATES THE INITIAL PARTICLES ★
  - Makes Deuterium particles with thermal velocities
  - Makes Tritium particles with thermal velocities
  - Positions them randomly inside the plasma boundary
  - Returns vector of all particles

**This is the file that generates your particles!**

### 5. **main.cpp** - Main Program
**What it does:**
- Sets up OpenGL window
- Creates tokamak geometry object
- Creates magnetic field object
- Creates plasma physics engine
- **CALLS `createThermalPlasma()` to get initial particles** (Line 113)
- Main loop:
  - Calls `updateParticles()` to evolve physics
  - Renders tokamak structure
  - Renders all particles
  - Prints fusion events to console

## WHERE AND HOW TO ADD PARTICLES

### OPTION 1: Modify Number of Initial Particles (EASIEST)

**Location:** `main.cpp`, lines 110-111

```cpp
// Change these numbers to create more or fewer particles
int numDeuterium = 100;  // ← Change this (try 200, 300, 500)
int numTritium = 100;    // ← Change this
```

**Effect:** More particles = more collisions = more fusion = more realistic but slower

### OPTION 2: Add Particles During Simulation (ADVANCED)

**Location:** `main.cpp`, inside the main loop (around line 145)

```cpp
// After the physics update, you can add new particles
plasmaPhysics.updateParticles(particles, deltaTime);

// ADD YOUR CODE HERE to inject new particles:
if (currentTime > lastInjectionTime + 1.0) {  // Every 1 second
    // Create a new deuterium particle at specific position
    Particle newD = createParticle(Particle::DEUTERIUM, 
                                   0.1f, 0.0f,        // position (x, y)
                                   1e5f, 0.0f);       // velocity (vx, vy)
    particles.push_back(newD);
    lastInjectionTime = currentTime;
}
```

**Effect:** Simulates continuous fueling (like real reactors)

### OPTION 3: Customize Particle Distribution

**Location:** `plasma_physics.h`, function `createThermalPlasma()` (line 327)

**Current code:**
```cpp
std::uniform_real_distribution<float> pos_dist(-0.2f, 0.2f);
```

**Modify to change initial distribution:**
```cpp
// Example: Start all particles near center
std::uniform_real_distribution<float> pos_dist(-0.1f, 0.1f);

// Example: Ring distribution (annulus)
float angle = dist(rng) * 2 * M_PI;
float r = 0.15f;  // Fixed radius
x = r * cos(angle);
y = r * sin(angle);
```

### OPTION 4: Add Different Particle Types

**Location:** `main.cpp`, after line 113

```cpp
std::vector<Particle> particles = plasmaPhysics.createThermalPlasma(numDeuterium, numTritium);

// ADD ELECTRONS for quasi-neutrality
for (int i = 0; i < 50; ++i) {
    float x = (rand() / (float)RAND_MAX) * 0.4f - 0.2f;
    float y = (rand() / (float)RAND_MAX) * 0.4f - 0.2f;
    if (tokamak.isInsidePlasma(x, y)) {
        Particle electron = createParticle(Particle::ELECTRON, x, y, 1e6f, 1e6f);
        particles.push_back(electron);
    }
}

// ADD IMPURITIES (e.g., Helium ash from previous fusions)
for (int i = 0; i < 10; ++i) {
    float x = (rand() / (float)RAND_MAX) * 0.4f - 0.2f;
    float y = (rand() / (float)RAND_MAX) * 0.4f - 0.2f;
    if (tokamak.isInsidePlasma(x, y)) {
        Particle helium = createParticle(Particle::HELIUM, x, y, 5e5f, 5e5f);
        particles.push_back(helium);
    }
}
```

## HOW THE PHYSICS FLOWS

```
1. main.cpp creates particles via plasmaPhysics.createThermalPlasma()
                    ↓
2. Each frame, main.cpp calls plasmaPhysics.updateParticles(particles, dt)
                    ↓
3. updateParticles() loops through all particles:
   - Calls applyMagneticForce() for each particle
       → This uses magnetic_field.getTotalField(x, y)
       → Calculates F = q(v × B)
       → Updates particle velocity
   
   - Calls applyCoulombForce() for each pair
       → Calculates F = k·q₁·q₂/r² with Debye screening
       → Updates both particles' velocities
   
   - Calls attemptFusion() for each D-T pair
       → Checks distance and energy
       → If conditions met, creates He + n products
       → Adds new particles to vector
   
   - Updates position: x += vx·dt, y += vy·dt
   
   - Calls checkBoundaryCollision()
       → Uses tokamak.distanceFromPlasmaEdge()
       → Reflects particles at wall
                    ↓
4. main.cpp renders all active particles
```

## SUMMARY: THE TWO MAIN PLACES TO TOUCH

### FOR INITIAL PARTICLES:
**File:** `main.cpp`
**Line:** 110-113
**What to change:** `numDeuterium` and `numTritium` values

### FOR PARTICLE BEHAVIOR:
**File:** `plasma_physics.h`
**Function:** `createThermalPlasma()`
**What to change:** Velocity distributions, positions, types

## EXAMPLE: ADDING A PARTICLE BEAM

Want to simulate a neutral beam injector? Add this to `main.cpp`:

```cpp
// After line 113 (after initial plasma creation)

// Inject high-energy deuterium beam from the side
for (int i = 0; i < 20; ++i) {
    float y_pos = -0.2f + i * 0.02f;  // Vertical spread
    Particle beamParticle = createParticle(
        Particle::DEUTERIUM,
        -0.3f, y_pos,      // Start from left side
        5e6f, 0.0f         // High horizontal velocity (5×10⁶ m/s)
    );
    particles.push_back(beamParticle);
}
```

## DEBUGGING TIPS

**Print particle info:**
```cpp
// In main.cpp, inside render loop
if (frameCount % 60 == 0) {  // Every 60 frames
    std::cout << "Active particles: " << particles.size() << std::endl;
    std::cout << "First particle pos: (" << particles[0].x << ", " 
              << particles[0].y << ")" << std::endl;
}
```

**Visualize forces:**
```cpp
// In plasma_physics.h, applyMagneticForce()
std::cout << "Particle at (" << p.x << ", " << p.y << ") "
          << "Force: (" << Fx << ", " << Fy << ")" << std::endl;
```

## PERFORMANCE NOTES

- **100 particles**: Real-time, smooth
- **200 particles**: Real-time, slight slowdown
- **500 particles**: Noticeable lag (O(N²) Coulomb interactions)
- **1000+ particles**: Very slow without optimization

**To handle more particles:**
1. Reduce Coulomb interaction range
2. Use spatial hashing (divide space into cells)
3. Implement GPU acceleration
4. Use particle-mesh methods for long-range forces

---

You now have everything you need! The particles are created in `plasma_physics.h::createThermalPlasma()` and updated in `plasma_physics.h::updateParticles()`. Modify the numbers in `main.cpp` to change how many you start with.
