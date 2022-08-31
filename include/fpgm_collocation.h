/*
* fpgm_collocation.h
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


// Flat Plate Glider Model with collocation functions

#ifndef FPGM_COLLOCATION_H
#define FPGM_COLLOCATION_H

#include <math.h>
#include <fstream>
#include <vector>
#include <string>

#include "math.hpp"
#include "geo.h"
#include "mathlib.h"
#include "Eigen/Dense"
#include <nlopt.h>

// https://stackoverflow.com/questions/5693686/how-to-use-yaml-cpp-in-a-c-program-on-linux
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace Eigen;
// using matrix::Vector2d;
// using matrix::Vector2f;
using matrix::wrap_pi;
using matrix::wrap_2pi;

namespace fpgm_collocation
{
    /** @brief Robust Post-Stall Perching with a Simple Fixed-Wing Glider using LQR-Trees
     * source : https://groups.csail.mit.edu/robotics-center/public_papers/Moore14a.pdf
     * orginal source : estimation of fpgm https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.208.1676&rep=rep1&type=pdf
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
     * x = [x, z, theta, phi, xdot, zdot, thetadot]
     * 7 dynamics,
     * dx = [xdot, zdot, thetadot, phidot, xdotdot, zdotdot, thetadotdot]
     * 1 input
     * u = [phidot]
     * 
    **/

    class equations_and_helper
    {
        public:
            struct fpgm_param
            {
                double l_w, l_e, l;
                double s_w, s_e; // Surface area of the wing and tail control surfaces
                double mass;
                double I; // I is only rotation in single axis
                double h; // Time-step
                Eigen::MatrixXd Q;
                double R;
            };

            struct optimization_constrain
            {
                double v_c;
                double t_c;
                double p_c;
                double td_c;
                double pd_c;
            };

            double cl(double aoa) { return 2 * sin(aoa) * cos(aoa);};
            
            double cd(double aoa) { return 2 * pow(sin(aoa), 2);};

            inline double two_d_cross(Eigen::Vector2d v1, Eigen::Vector2d v2) {return v1.x() * v2.y() - v1.y() * v2.x();};

            /** @brief calculate the dynamics (time_derivative of state)
             * @param x 
             * @param z 
             * @param theta 
             * @param phi
             * @param xdot
             * @param zdot
             * @param thetadot
             * @param phidot
             * **/
            Eigen::VectorXd fpgm_dynamics(
                double x, double z, double theta, double phi, double xdot, double zdot, double thetadot, double phidot,
                fpgm_param parameter)
            {
                double g = 9.81 , p = 1.225; // Density of air = 1.225 kg/m

                Eigen::VectorXd dx(7);

                // force_vectors
                Eigen::Vector2d n_w = Eigen::Vector2d(-sin(theta), cos(theta));
                Eigen::Vector2d n_e = Eigen::Vector2d(-sin(theta + phi), cos(theta + phi));

                Eigen::Vector2d x_w = Eigen::Vector2d(
                    x - parameter.l_w * cos(theta), 
                    z - parameter.l_w * sin(theta));
                Eigen::Vector2d x_e = Eigen::Vector2d(
                    x - parameter.l * cos(theta) - parameter.l_e * cos(theta + phi), 
                    z - parameter.l * sin(theta) - parameter.l_e * sin(theta + phi));
                
                Eigen::Vector2d x_w_dot = Eigen::Vector2d(
                    xdot + parameter.l_w * thetadot * sin(theta), 
                    zdot - parameter.l_w * thetadot * cos(theta));
                Eigen::Vector2d x_e_dot = Eigen::Vector2d(
                    xdot + parameter.l * thetadot * sin(theta) + parameter.l_e * (thetadot + phidot) * sin(theta + phi), 
                    zdot - parameter.l * thetadot * cos(theta) - parameter.l_e * (thetadot + phidot) * cos(theta + phi));

                double alpha_w = theta - atan(x_w_dot[1] / x_w_dot[0]);
                double alpha_e = theta + phi - atan(x_e_dot[1] / x_e_dot[0]);
                Eigen::Vector2d force_w = 
                    0.5 * p * pow(pow(x_w_dot[0],2) + pow(x_w_dot[1],2), 2) * 
                    parameter.s_w * (cl(alpha_w) + cd(alpha_w)) * n_w;
                Eigen::Vector2d force_e =
                    0.5 * p * pow(pow(x_e_dot[0],2) + pow(x_e_dot[1],2), 2) * 
                    parameter.s_e * (cl(alpha_e) + cd(alpha_e)) * n_e;

                Eigen::Vector2d pos_dotdot = (force_w + force_e - Eigen::Vector2d(0, parameter.mass * g)) / parameter.mass;
                double theta_dotdot = 
                    (two_d_cross(Eigen::Vector2d(parameter.l_w, 0), (force_w)) + 
                    two_d_cross(Eigen::Vector2d(-parameter.l - parameter.l_e * cos(theta), 
                    -parameter.l + parameter.l_e * sin(theta)), (force_e))) / parameter.I;

                // dx = [xdot, zdot, thetadot, phidot, xdotdot, zdotdot, thetadotdot]
                dx[0] = xdot;
                dx[1] = zdot;
                dx[2] = thetadot;
                dx[3] = phidot;
                dx[4] = pos_dotdot[0];
                dx[5] = pos_dotdot[1];
                dx[6] = theta_dotdot;

                return dx;
            }

            Eigen::VectorXd std_vector_to_eigen_vector(std::vector<double> x)
            {
                int vector_size = (int)x.size();
                Eigen::VectorXd v(vector_size);
                for (int i = 0; i < vector_size; i++)
                    v[i] = x[i];
                return v;
            }

    };

    class fpgm_collocation
    {

        private:
            
            equations_and_helper::fpgm_param param;
            equations_and_helper::optimization_constrain boundary;
            int N;

            std::vector<double> guess;

            /** @brief Collocation method used is the trapezoidal collocation
             * reference : https://epubs.siam.org/doi/pdf/10.1137/16M1062569
             * 
             * @param x = vector of all the states compressed into 1 dimension
             * Total size is 8 variables * N steps
             *  
            **/
            static void collocation_eq_constraints(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data)
            {
                // Data will provide start state

                static equations_and_helper eq;
                equations_and_helper::fpgm_param *params = (equations_and_helper::fpgm_param*)data;

                int equality_size = m / 7;

                for (int i = 0; i < equality_size; i++)
                {
                    // current state 
                    std::vector<double> x1 {
                        x[0+8*i], x[1+8*i], x[2+8*i], x[3+8*i], 
                        x[4+8*i], x[5+8*i], x[6+8*i]};
                    // future state
                    std::vector<double> x2 {
                        x[0+8*(i+1)], x[1+8*(i+1)], x[2+8*(i+1)], x[3+8*(i+1)], 
                        x[4+8*(i+1)], x[5+8*(i+1)], x[6+8*(i+1)]};

                    // current dynamics 
                    Eigen::VectorXd f_k = eq.fpgm_dynamics(
                        x1[0], x1[1], x1[2], x1[3], 
                        x1[4], x1[5], x1[6], x[7+8*i],
                        *params);
                    // future dynamics
                    Eigen::VectorXd f_k_1 = eq.fpgm_dynamics(
                        x2[0], x2[1], x2[2], x2[3], 
                        x2[4], x2[5], x2[6], x[7+8*(i+1)],
                        *params);

                    Eigen::VectorXd x_k = eq.std_vector_to_eigen_vector(x1);
                    Eigen::VectorXd x_k_1 = eq.std_vector_to_eigen_vector(x2);
                    
                    // Eigen::VectorXd lhs = 0.5 * params->h * (f_k_1 + f_k);
                    Eigen::VectorXd lhs = params->h * f_k;
                    Eigen::VectorXd rhs = x_k_1 - x_k;

                    Eigen::VectorXd single_results_vector = lhs - rhs;

                    for (int j = 0; j < 7; j++)
                        result[j+i*7] = single_results_vector[j];
                }
            }


            static double control_effort_objective(unsigned n, const double *x, double *grad, void *data)
            {
                equations_and_helper::fpgm_param *params = (equations_and_helper::fpgm_param*)data;

                // Assuming h_k = uniform, timestep is uniform
                // double factor = params->h / 2;
                double factor = params->h;
                double cost = 0;
                int state_input_length = n / 8;

                // for (int i = 0; i < state_input_length; i++)
                //     cost += factor * (pow(x[7+8*i], 2) + pow(x[7+8*(i+1)], 2));

                for (int i = 0; i < state_input_length; i++)
                {
                    Eigen::VectorXd x1(7);
                    x1 << x[0+8*i], x[1+8*i], x[2+8*i], x[3+8*i], 
                        x[4+8*i], x[5+8*i], x[6+8*i];
                    
                    auto state_term = x1.transpose() * params->Q * x1;
                    
                    // cout << x1.transpose() * params->Q << endl;
                    // cout << x1.transpose() << endl;
                    // cout << state_term << endl;

                    double input_term = x[7+8*i] * params->R * x[7+8*i];
                    cost += state_term + input_term;
                }

                cost *= factor;

                printf("cost = %lf\n", cost);
                return cost;
            }

        public:
            

            bool load_parameters(std::string directory, double total, int size, MatrixXd Q, double R)
            {
                ifstream f(directory.c_str());
                if (!f.good())
                    return false;

                YAML::Node node = YAML::LoadFile(directory);

                param = {}; boundary = {}; // reset the parameters
                param.l_w = node["length_cg_to_cwing"].as<double>();
                param.l_e = node["length_pivote_to_celevator"].as<double>();
                param.l = node["length_cg_to_pivote"].as<double>();
                param.s_e = node["surface_area_elevator"].as<double>();
                param.s_w = node["surface_area_wing"].as<double>();
                param.mass = node["mass"].as<double>();
                param.I = node["moments_of_inertia"].as<double>();
                param.Q = Q;
                param.R = R;
                param.h = total / (size-1);

                boundary.v_c = node["velocity_constrain"].as<double>();
                boundary.t_c = node["theta_contrain"].as<double>();
                boundary.p_c = node["phi_contrain"].as<double>();
                boundary.td_c = node["thetadot_constrain"].as<double>();
                boundary.pd_c = node["phidot_constrain"].as<double>();

                printf("Parameters loaded\n");
                return true;
            }

            bool load_initial_guess(std::vector<double> x)
            {
                guess.clear();
                guess = x;
                if (remainder((int)x.size(), 8))
                    return false;
                N = ((int)x.size() / 8);
                printf("guess size = %d, N steps = %d\n", (int)guess.size(), N);
                return true;
            }

            bool set_bounds(double *lb, double *ub)
            {
                int size_of_vector = 
                    (sizeof(lb) / sizeof(double)) / 8;
                for (int i = 0; i < size_of_vector; i++)
                {
                    // Give a threshold that is from the guess

                    lb[0+i*8] = -1E8; ub[0+i*8] = 1E8; // x constrain
                    lb[1+i*8] = -1E8; ub[1+i*8] = 1E8; // z constrain
                    
                    lb[2+i*8] = -boundary.t_c; ub[2+i*8] = boundary.t_c; // theta constrain
                    lb[3+i*8] = -boundary.p_c; ub[3+i*8] = boundary.p_c; // phi constrain
                    
                    // velocity x and z constrain
                    lb[4+i*8] = lb[5+i*8] = -boundary.v_c;
                    ub[4+i*8] = ub[5+i*8] = boundary.v_c;

                    lb[6+i*8] = -boundary.td_c; ub[6+i*8] = boundary.td_c; // thetadot constrain
                    lb[7+i*8] = -boundary.pd_c; ub[7+i*8] = boundary.pd_c; // phidot constrain
                }

                return true;
            }

            Eigen::Vector3d differential_flat_estimated_rotation(Eigen::Vector3d a, double y)
            {
                Eigen::Vector3d alpha = a + Eigen::Vector3d(0,0,9.81);
                Eigen::Vector3d xC(cos(y), sin(y), 0);
                Eigen::Vector3d yC(-sin(y), cos(y), 0);
                Eigen::Vector3d xB = (yC.cross(alpha)).normalized();
                Eigen::Vector3d yB = (alpha.cross(xB)).normalized();
                Eigen::Vector3d zB = xB.cross(yB);

                Eigen::Matrix3d R;
                R.col(0) = xB;
                R.col(1) = yB;
                R.col(2) = zB;

                double yaw = atan2(R(1,0), R(0,0));
                double pitch = atan2(-R(2,0), sqrt(pow(R(2,1),2) + pow(R(2,2),2)));
                double roll = atan2(R(2,1), R(2,2));

                return Eigen::Vector3d(roll, pitch, yaw);
            }

            Eigen::VectorXd nlopt_optimization() 
            {
                Eigen::VectorXd final_vector;
                if (guess.empty())
                    return final_vector;
                
                int equality_vector_size = N*(8-1);
                double tol[equality_vector_size];

                nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, guess.size());
                nlopt_set_min_objective(opt, control_effort_objective, &param);
                
                nlopt_set_ftol_abs(opt, 1E-6);
                nlopt_set_xtol_rel(opt, 1E-3);
                nlopt_set_maxeval(opt, 1E2);
                nlopt_set_maxtime(opt, 2.0); 

                // nlopt_add_equality_mconstraint(
                //     opt, equality_vector_size, collocation_eq_constraints, &param, NULL);
                nlopt_add_equality_mconstraint(
                    opt, equality_vector_size, collocation_eq_constraints, &param, tol);

                double lb[guess.size()], ub[guess.size()];
                set_bounds(lb, ub);

                // nlopt_set_lower_bounds(opt, lb);
                // nlopt_set_upper_bounds(opt, ub);

                double x[guess.size()];
                std::copy(guess.begin(), guess.end(), x);

                int x_size = sizeof(x) / sizeof(int);
                // printf("guess_size = %d\n", (int)guess.size());
                // printf("x:\n");
                // for (int i = 0; i < guess.size(); i++)
                //     printf("%lf ", x[i]);
                // printf("\n");

                double cost = 0;
                nlopt_optimize(opt, x, &cost);
                printf("number of iterations: %d \n", nlopt_get_numevals(opt));

                printf("guess-difference: \n");
                for(int i = 0; i < N; i++)
                {
                    printf("row %d ", i);
                    for (int j = 0; j < 8; j++)
                        printf("%lf ", x[j+i*8] - guess[j+i*8]);
                    printf("\n");
                }
                printf("\n");

                printf("Optimization completed cost %lf\n", cost);

            }
    };

}

#endif