/*
* opt_landing.cpp
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
#include <fstream>
#include <string>
#include <chrono>
#include <vector>

#include "obvp.h"
#include "fpgm_collocation.h"
#include "matplotlibcpp.h"

// https://stackoverflow.com/questions/5693686/how-to-use-yaml-cpp-in-a-c-program-on-linux
#include "yaml-cpp/yaml.h"

using namespace obvp;
using namespace fpgm_collocation;
using namespace std::chrono;
namespace plt = matplotlibcpp;

double rad_to_deg = 1/3.14159265358979323846264338327 * 180.0;
double deg_to_rad = 1/180.0 * 3.14159265358979323846264338327;
matrix::Vector3d zero = matrix::Vector3d{0, 0, 0};

// Don't comprehend this 
// since it is in ENU frame, hence it is pry not rpy
matrix::Vector3d velocity_in_global_frame(double y, double r, double p, double vel) 
{
    // in the enu frame, x would be the velocity
    matrix::Vector3d relative_vel = matrix::Vector3d(vel, 0, 0);
    double yaw[9] = {
            1.0,  0.0,  0.0,
            0.0,  cos(y),  -sin(y),
            0.0,  sin(y),  cos(y)
        };
    double roll[9] = {
            cos(r),  0.0,  sin(r),
            0.0,  1.0,  0.0,
            -sin(r),  0.0,  cos(r)
        };
    
    double pitch[9] = {
            cos(p),  -sin(p), 0.0,
            sin(p),  cos(p), 0.0,
            0.0,  0.0,  1.0
        };
    matrix::SquareMatrix<double, 3> P(pitch);
    matrix::SquareMatrix<double, 3> R(roll);
    matrix::SquareMatrix<double, 3> Y(yaw);
    return (P * R * Y) * relative_vel; 
}



int main(int argc, char **argv) 
{
    MapProjection _global_local_proj_ref{};

    fpgm_collocation::fpgm_collocation fpgm;

    std::string params_directory = "parameters.yaml";

    ifstream f(params_directory.c_str());
    if (!f.good())
        return false;

    YAML::Node node = YAML::LoadFile(params_directory);
    
    /**
     * @brief Parameters setting
     */
    // ENU frame according to Lat Lon and height
    double WS_LAT_P_LAND = 1.330587;
    double WS_LONG_P_LAND = 103.783740;

    // Setting the reference of the map with the current position
    _global_local_proj_ref.initReference(WS_LAT_P_LAND, WS_LONG_P_LAND);

    double airspeed = node["airspeed"].as<double>(); // m/s
    double descend_bearing_deg = node["descend_bearing_deg"].as<double>(); // bearing in deg
    double descend_pitch_deg = node["descend_pitch_deg"].as<double>(); // angle in deg
    double buffer_distance = node["buffer_distance"].as<double>(); // angle in deg
    double descend_bearing_rad = descend_bearing_deg * deg_to_rad;
    double descend_pitch_rad = descend_pitch_deg * deg_to_rad;

    double height_of_descend = node["height_of_descend"].as<double>(); // m
    double height_of_land = node["height_of_land"].as<double>(); // If there is a structure for the aircraft to land on

    // according to Robust Post-Stall Perching with a Fixed-Wing UAV by Joseph Moore page 28
    // elevator moving upwards is considered negative
    double max_elevator = node["phi_contrain"].as<double>();
    double max_elevator_rad = max_elevator * deg_to_rad;
    double current_elevator_rad = node["current_elevator_rad"].as<double>();
    double current_thetadot = node["current_thetadot"].as<double>();

    double command_time = node["command_time"].as<double>();

    // Calculate the descend distance
    double distance_to_land_from_dive = 
        (height_of_descend - height_of_land) / tan(descend_pitch_rad) + buffer_distance;
    printf("distance_to_land_from_dive = %lf\n", distance_to_land_from_dive);

    matrix::Vector3d dive_position;
    matrix::Vector3d landing_position = matrix::Vector3d(
        WS_LAT_P_LAND, WS_LONG_P_LAND, height_of_land);
    matrix::Vector3d platform_speed;
    dive_position(2) = height_of_descend;

    double descend_bearing_backwards = wrap_pi(descend_bearing_rad-3.14);
    waypoint_from_heading_and_distance(
        WS_LAT_P_LAND, WS_LONG_P_LAND, descend_bearing_backwards, 
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


    /**
     * @brief require estimates for x(state) = [x, z, theta, phi, xdot, zdot, thetadot]
     * @arg x, z, xdot, zdot can be taken care of by the bvp solution
     * @arg thetadot can be taken care of, and we assume it to be constant linear interpolation from <current point> to 0
     * @arg phi can be taken care of, and we assume it to be constant linear interpolation from <current point> to max_phi
     * @arg theta can be assumed with differentially flat outputs from the bvp
     */

    /** @brief Boundary value problem provides x, z, xdot, zdot estimates */
    matrix::Vector3d alpha, beta, gamma;

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

    // Make a guess to reach the final goal point
    // distance / sqrt(vx^2 + vy^2) * factor
    double guess_factor = 3.0;
    double stepping_factor = 1.0/10.0;

    // make sure total time is a factor of command time
    double total_time_division = 
        round(distance_to_land_from_dive / 
        sqrt(pow(velocity_global(0),2) + pow(velocity_global(1),2)) / 
        command_time);

    double total_time = total_time_division * command_time * guess_factor;
    // BVP guessing step (to reduce by)
    double step = command_time * guess_factor;
    bool check_passed = false;
    int iter = 0;
    px4_array_container waypoints;

    time_point<std::chrono::system_clock> bvp_start = system_clock::now();
    while (!check_passed)
    {
        iter += 1;
        get_bvp_coefficients(initial_state_local, final_state_local, total_time,
            &alpha, &beta, &gamma);

        int bad_counts = check_z_vel(
            initial_state_local, final_state_local, total_time, command_time, 
            alpha, beta, gamma);
        
        if (bad_counts == 0)
            check_passed = true;
        else
        {
            // printf("iter %d with bad_counts : %d\n", iter, bad_counts);
            total_time -= (double)bad_counts * stepping_factor * step;
        }
    }
    
    auto bvp_time = duration<double>(system_clock::now() - bvp_start).count();
    printf("bvp_time taken : %lfs with total calc time : %lfs\n", bvp_time, total_time);

    int waypoint_size = 0;
    waypoints = get_discrete_points(
        initial_state_local, final_state_local, total_time, command_time, 
        alpha, beta, gamma, waypoint_size);

    // printf("discrete waypoints with size %d :\n", waypoint_size);
    // for (int i = 0; i < waypoint_size; i++)
    //     printf("[%lf %lf %lf]\n", waypoints[i](0), waypoints[i](1), waypoints[i](2));
    // printf("\n");

    // printf("alpha coefficient(%lf, %lf, %lf)\n", alpha(0), alpha(1), alpha(2));
    // printf("beta coefficient(%lf, %lf, %lf)\n", beta(0), beta(1), beta(2));
    // printf("gamma coefficient(%lf, %lf, %lf)\n", gamma(0), gamma(1), gamma(2));
    // printf("waypoint_size(%d) iter(%d)\n", (int)waypoints.size(), iter);
    
    /** @brief theta estimates */
    vector<double> theta_vector, phi_vector, thetadot_vector; // pitch
    double phi_factor = max_elevator_rad / (double)(waypoint_size-1);

    double trim = - descend_pitch_rad;
    for (int i = 0; i < waypoint_size; i++)
    {
        theta_vector.push_back(
            trim - fpgm.differential_flat_estimated_rotation(
            Eigen::Vector3d(waypoints[i](6), waypoints[i](7), waypoints[i](8)), 
            descend_bearing_rad)[1]);
        phi_vector.push_back(current_elevator_rad - phi_factor*i);
        thetadot_vector.push_back(current_thetadot);
    }

    /** @brief do a transformation to orientate waypoints to align with x axis */
    std::vector<matrix::Vector2d> vector_t_waypoints, vector_t_velocity;
    double yaw[4] = {
        cos(descend_bearing_backwards),  -sin(descend_bearing_backwards),
        sin(descend_bearing_backwards),  cos(descend_bearing_backwards)
    };
    matrix::SquareMatrix<double, 2> Y(yaw);
    for (int i = 0; i < waypoint_size; i++)
    {
        matrix::Vector2d t_waypoints = Y.I() *
            matrix::Vector2d(waypoints[i](0), waypoints[i](1));
        matrix::Vector2d t_velocity = Y.I() *
            matrix::Vector2d(waypoints[i](3), waypoints[i](4));
        
        vector_t_waypoints.push_back(t_waypoints);
        vector_t_velocity.push_back(t_velocity);
        printf("pos [%lf %lf] vel [%lf %lf] theta [%lf] phi [%lf]\n", 
            t_waypoints(0), waypoints[i](2),
            t_velocity(0), waypoints[i](5),
            theta_vector[i], phi_vector[i]);
    }

    std::vector<double> initial_guess;
    std::vector<double> initial_x, initial_z;
    for (int i = 0; i < waypoint_size; i++)
    {
        // x = [x, z, theta, phi, xdot, zdot, thetadot]
        // u = [phidot]
        initial_guess.push_back(vector_t_waypoints[i](0));
        initial_guess.push_back(waypoints[i](2));
        initial_guess.push_back(theta_vector[i]);
        initial_guess.push_back(phi_vector[i]);
        initial_guess.push_back(vector_t_velocity[i](0));
        initial_guess.push_back(waypoints[i](5));
        initial_guess.push_back(thetadot_vector[i]);
        initial_guess.push_back(0.0);

        initial_x.push_back(vector_t_waypoints[i](0));
        initial_z.push_back(waypoints[i](2));
    }
    Eigen::Matrix< double, 7, 1> v;
    v << 0.05, 0.10, 300.0, 300.0, 0.05, 0.05, 2000.0;
    Eigen::Matrix< double, 7, 7> Q = v.array().matrix().asDiagonal();
    // cout << Q << endl;
    double R = 200;

    if (!fpgm.load_parameters(
        params_directory, total_time, 
        (int)initial_guess.size(), Q, R,
        initial_x, initial_z))
        return -1;

    if (!fpgm.load_initial_guess(initial_guess))
        return -1;
    
    time_point<std::chrono::system_clock> opt_start = system_clock::now();
    fpgm_collocation::fpgm_collocation::control_state control;
    control = fpgm.nlopt_optimization();
    auto opt_time= duration<double>(system_clock::now() - opt_start).count();
    printf("opt_time taken : %lfs\n", opt_time);


    /** @brief Visualization **/
    // Set the size of output image to 1200x780 pixels
    plt::figure_size(980, 460);
    // plot a red dashed line from given x and y data.
    plt::named_plot("optimal", control.x, control.z, "r--");
    plt::named_plot("guess", initial_x, initial_z, "b--");
    
    // Enable legend.
    plt::legend();
    
    string title = "height against x";
    plt::title(title); // add graph title
    plt::grid(true);
    plt::show();

    return 0;
}