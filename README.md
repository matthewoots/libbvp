# Optimal Boundary Value Problem Library(libbvp)

## Installation
`libbvp` optimal boundary value problem library using some PX4 modules. OBVP is referenced from Mark Mueller from https://flyingmachinearena.org/wp-content/uploads/mueTRO15.pdf

### PX4 Libraries used
- mathlib : Some useful definitions for math classes
- PX4-Matrix : To use the vector classes that are similar to Eigen library, even multiplication is supported (more info at https://github.com/PX4/PX4-Autopilot/tree/4a3d64f1d76856d22323d1061ac6e560efda0a05/src/lib/matrix)
- geo : To use gps coordinates (lat and lon) to identify local coordinates or targets (https://github.com/PX4/PX4-Autopilot/tree/master/src/lib/geo)

### Setup
```bash
git clone git@github.com:matthewoots/libbvp.git --recurse-submodules
mkdir build && cd build
cmake ..
make
```