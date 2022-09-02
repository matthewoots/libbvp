/*
* bvp_test.cpp
*
* ---------------------------------------------------------------------
* Copyright (C) 2022 Matthew (matthewoots at gmail.com)
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License
*  as published by the Free Software Foundation; either version 2
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
* ---------------------------------------------------------------------
*/

#include <iostream>
#include <string>
#include "obvp.h"

using namespace obvp;

double rad_to_deg = 1/3.14159265358979323846264338327 * 180.0;
double deg_to_rad = 1/180.0 * 3.14159265358979323846264338327;
matrix::Vector3d zero = matrix::Vector3d{0, 0, 0};

// since it is in ENU frame, hence it is pry not rpy
matrix::Vector3d velocity_in_global_frame(double r, double p, double y, double vel) 
{
    // in the enu frame, y would be the velocity
    matrix::Vector3d relative_vel = matrix::Vector3d(0, vel, 0);
    double pitch[9] = {
            1.0,  0.0,  0.0,
            0.0,  cos(p),  -sin(p),
            0.0,  sin(p),  cos(p)
        };
    matrix::SquareMatrix<double, 3> P(pitch);
    double yaw[9] = {
            cos(y),  -sin(y), 0.0,
            sin(y),  cos(y), 0.0,
            0.0,  0.0,  1.0
        };
    matrix::SquareMatrix<double, 3> Y(yaw);
    return (Y.I() * P.I()) * relative_vel; 
}

/**
 * @brief bvp test case
 */

int main(int argc, char **argv) 
{
    MapProjection _global_local_proj_ref{};
    
    // ENU frame according to Lat Lon and height
    double WS_LAT_P_LAND = 1.330587;
    double WS_LONG_P_LAND = 103.783740;
    _global_local_proj_ref.initReference(WS_LAT_P_LAND, WS_LONG_P_LAND);

    double airspeed = 20.0; // m/s
    double descend_bearing_deg = 10.0; // bearing in deg
    double descend_pitch_deg = 30.0; // angle in deg
    double descend_bearing_rad = descend_bearing_deg * deg_to_rad;
    double descend_pitch_rad = descend_pitch_deg * deg_to_rad;

    double height_of_descend = 10.0; // m
    double height_of_land = 0.5; // If there is a structure for the aircraft to land on

    double distance_to_land_from_dive = (height_of_descend - height_of_land) / sin(descend_pitch_rad);

    matrix::Vector3d dive_position;
    matrix::Vector3d landing_position = matrix::Vector3d(
        WS_LAT_P_LAND, WS_LONG_P_LAND, height_of_land);
    matrix::Vector3d platform_speed;
    dive_position(2) = height_of_descend;

    waypoint_from_heading_and_distance(
        WS_LAT_P_LAND, WS_LONG_P_LAND, -descend_bearing_rad, 
        (float)distance_to_land_from_dive, &dive_position(0), &dive_position(1));

    matrix::Vector3d velocity_global = 
        velocity_in_global_frame(0.0, descend_pitch_rad, descend_bearing_rad, airspeed);

    matrix::Vector3d dive_position_local;
    matrix::Vector3d landing_position_local = 
        matrix::Vector3d(0, 0, landing_position(2));
    if (_global_local_proj_ref.isInitialized()) 
    {
        matrix::Vector2f pos_local = _global_local_proj_ref.project(dive_position(0), dive_position(1));
        dive_position_local(0) = (double)pos_local(0);
        dive_position_local(1) = (double)pos_local(1);
        dive_position_local(2) = dive_position(2);
    }

    printf("dive_position [%lf %lf %lf] dive_position_local [%lf %lf %lf]\n", 
        dive_position(0), dive_position(1), dive_position(2),
        dive_position_local(0), dive_position_local(1), dive_position_local(2));
    printf("landing_position [%lf %lf %lf] landing_position_local [%lf %lf %lf]\n", 
        landing_position(0), landing_position(1), landing_position(2),
        landing_position_local(0), landing_position_local(1), landing_position_local(2));
    printf("velocity [%lf %lf %lf]\n", velocity_global(0), velocity_global(1), velocity_global(2));

    matrix::Vector3d alpha, beta, gamma;
    double total_time = 5.0;
    matrix::SquareMatrix<double, 3> initial_state_local;
    initial_state_local.setZero();
    
    initial_state_local.col(0) = dive_position_local; // Position
    initial_state_local.col(1) = velocity_global; // Velocity
    initial_state_local.col(2) = zero; // Acceleration

    matrix::SquareMatrix<double, 3> final_state_local;
    final_state_local.setZero();
    
    final_state_local.col(0) = landing_position_local; // Position
    final_state_local.col(1) = zero; // Velocity
    final_state_local.col(2) = zero; // Acceleration

    get_bvp_coefficients(initial_state_local, final_state_local, total_time,
		&alpha, &beta, &gamma);

    printf("alpha coefficient(%lf, %lf, %lf)\n", alpha(0), alpha(1), alpha(2));
    printf("beta coefficient(%lf, %lf, %lf)\n", beta(0), beta(1), beta(2));
    printf("gamma coefficient(%lf, %lf, %lf)\n", gamma(0), gamma(1), gamma(2));

    return 0;
}