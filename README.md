# Precision Landing Optimal Boundary Value Problem Library(libprecisionbvp)

## Installation
`libprecisionbvp` optimal boundary value problem library using some PX4 modules. OBVP is referenced from Mark Mueller from https://flyingmachinearena.org/wp-content/uploads/mueTRO15.pdf

In addition `fpgm_collocation` header file consist of the flat plane glider dynamics and non-linear optimization solving the optimal trajectory problem for a `fpgm` type of fixed wing.

This precision landing maneuver utilizes the Boundary Value Problem to provide a guess for the optimization and direct collocation method used in propagating the states forward in time.

### PX4 Libraries used
- mathlib : Some useful definitions for math classes
- PX4-Matrix : To use the vector classes that are similar to Eigen library, even multiplication is supported (more info at https://github.com/PX4/PX4-Autopilot/tree/4a3d64f1d76856d22323d1061ac6e560efda0a05/src/lib/matrix)
- geo : To use gps coordinates (lat and lon) to identify local coordinates or targets (https://github.com/PX4/PX4-Autopilot/tree/master/src/lib/geo)
- eigen : From PX4 repository (https://github.com/PX4/eigen) - seems like eigen could be used all along

## Dependencies
- Using `Nlopt` for non-linear optimization
```bash
cd
git clone git@github.com:stevengj/nlopt.git
cd nlopt
# commit that I used in this repository
git reset --hard 401330a359d53f1794ea211f23b4a272234b4787 
mkdir build
cd build
cmake ..
make
sudo make install
```

### Setup
```bash
git clone git@github.com:matthewoots/libbvp.git --recurse-submodules
mkdir build && cd build
cmake ..
make
```
Run with `./obvp_precision_landing` or `./obvp_opt_landing`