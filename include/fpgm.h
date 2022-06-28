// Flat Plate Glider Model

#ifndef FPGM_H
#define FPGM_H

#include <math.h>
#include "math.hpp"
#include "geo.h"
#include "mathlib.h"
#include "Eigen/Dense"

using namespace std;
using matrix::Vector2d;
using matrix::Vector2f;
using matrix::wrap_pi;
using matrix::wrap_2pi;

namespace fpgm
{
    /** @brief Robust Post-Stall Perching with a Simple Fixed-Wing Glider using LQR-Trees
     * source : https://groups.csail.mit.edu/robotics-center/public_papers/Moore14a.pdf
     * @param (double)alpha = angle of attack (rads)
     * @param (double)phi = elevator control angle (rads)
     * @param (double)theta = pitch angle (rads)
     * @param (matrix::Vector2d)Fw = force on wing
     * @param (matrix::Vector2d)Fe = force on elevator
     * @param (matrix::Vector2d)Xw = kinematics of the geometric centroid for wing
     * @param (matrix::Vector2d)Xe = kinematics of the geometric centroid for elevator
     * @param (matrix::Vector2d)nw = unit vectors normal to the wing control surface
     * @param (matrix::Vector2d)ne = unit vectors normal to the elevator control surface
     * 
     * 7 states, 
     * x = [x, z, θ , φ ,  ̇x,  ̇z,  ̇θ ]
     * 
     * **/

}

#endif