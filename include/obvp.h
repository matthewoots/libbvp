/*
* obvp.h
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

#ifndef OBVP_H
#define OBVP_H

#include <math.h>
#include <vector>
#include "math.hpp"
#include "geo.h"
#include "mathlib.h"
#include "Array.hpp"
#include "Eigen/Dense"

using namespace std;
using matrix::Vector2d;
using matrix::Vector2f;
using matrix::wrap_pi;
using matrix::wrap_2pi;
using px4_array_container = px4::Array< matrix::Vector<double, 9>,100>;
using pva = vector<matrix::Vector<double, 9>>;
using namespace Eigen;

namespace obvp
{
	// get_bvp_coefficients using PX4 matrix (without eigen)
    void get_bvp_coefficients(matrix::SquareMatrix<double, 3> initial,
		matrix::SquareMatrix<double, 3> final, double total_time,
		matrix::Vector3d *alpha, matrix::Vector3d *beta,
		matrix::Vector3d *gamma)
	{
		double T = total_time;
		// Update the initial values
		matrix::Vector3d p0 = matrix::Vector3d(
			initial(0,0), initial(1,0), initial(2,0));
		matrix::Vector3d v0 = matrix::Vector3d(
			initial(0,1), initial(1,1), initial(2,1));
		matrix::Vector3d a0 = matrix::Vector3d(
			initial(0,2), initial(1,2), initial(2,2));

		// Update the destination values
		matrix::Vector3d pf = matrix::Vector3d(
			final(0,0), final(1,0), final(2,0));
		matrix::Vector3d vf = matrix::Vector3d(
			final(0,1), final(1,1), final(2,1));
		matrix::Vector3d af = matrix::Vector3d(
			final(0,2), final(1,2), final(2,2));

		matrix::Vector<double, 9> delta;
		delta(0) = pf(0) - p0(0) - v0(0) * T - 0.5 * a0(0) * pow(T,2);
		delta(1) = pf(1) - p0(1) - v0(1) * T - 0.5 * a0(1) * pow(T,2);
		delta(2) = pf(2) - p0(2) - v0(2) * T - 0.5 * a0(2) * pow(T,2);
		delta(3) = (vf(0) - v0(0) - a0(0) * T);
		delta(4) = (vf(1) - v0(1) - a0(1) * T);
		delta(5) = (vf(2) - v0(2) - a0(2) * T);
		delta(6) = (af(0) - a0(0));
		delta(7) = (af(1) - a0(1));
		delta(8) = (af(2) - a0(2));

		double fifth_order[9] = {720, -360*T, 60*pow(T,2),
					-360*T, 168*pow(T,2), -24*pow(T,3),
					60*pow(T,2), -24*pow(T,3), 3*pow(T,4)};

		matrix::SquareMatrix<double, 3> m(fifth_order);

        // Print out the delta vector
        // printf("delta:\n");
        // for (int i = 0; i < 9; i++)
        // {
        //     printf("%lf ", delta(i));
        // }
        // printf("\n");

		matrix::SquareMatrix<double, 9> M;
		M.setZero();

		// Make the 3x3 matrix into a 9x9 matrix so that xyz is included
		for (int i = 0; i < 9; i++)
		{
			int l_m = (i+1) % 3 + 3*(int)(((i+1) % 3) == 0) - 1;
			int h_m = (int)ceil((double)(i+1) / 3.0) - 1;

			int l_M = (l_m)*3;
			int h_M = (h_m)*3;

			M(l_M,h_M) = m(l_m,h_m);
			M(l_M+1,h_M+1) = m(l_m,h_m);
			M(l_M+2,h_M+2) = m(l_m,h_m);

            // printf("M(%d,%d) to M(%d,%d)\n", l_M, l_M, l_M+2, h_M+2);
            // printf("m(%d,%d), %lf\n", l_m, h_m, m(l_m,h_m));
		}

        // Print out the M matrix vector
        // printf("M matrix:\n");
        // for (int i = 0; i < 9; i++)
        // {
        //     for (int j = 0; j < 9; j++)
        //     {
        //         printf("%lf ", M(i,j));
        //     }
        //     printf("\n");
        // }
        // printf("\n");

		matrix::Vector<double, 9> abg;

		abg = (1/pow(T,5) * M * delta);
		// printf("abg:\n");
        // for (int i = 0; i < 9; i++)
        //     printf("%lf ", abg(i));
        // printf("\n");

        *alpha = matrix::Vector3d(abg(0), abg(1), abg(2));
        *beta = matrix::Vector3d(abg(3), abg(4), abg(5));
        *gamma = matrix::Vector3d(abg(6), abg(7), abg(8));

	}

	// get_bvp_coefficients using eigen
	void get_bvp_coefficients(Eigen::Matrix3d initial,
		Eigen::Matrix3d final, double total_time,
		Eigen::Vector3d *alpha, Eigen::Vector3d *beta,
		Eigen::Vector3d *gamma)
	{
		double T = total_time;
		// Update the initial values
		Eigen::Vector3d p0 = Eigen::Vector3d(
			initial(0,0), initial(1,0), initial(2,0));
		Eigen::Vector3d v0 = Eigen::Vector3d(
			initial(0,1), initial(1,1), initial(2,1));
		Eigen::Vector3d a0 = Eigen::Vector3d(
			initial(0,2), initial(1,2), initial(2,2));

		// Update the destination values
		Eigen::Vector3d pf = Eigen::Vector3d(
			final(0,0), final(1,0), final(2,0));
		Eigen::Vector3d vf = Eigen::Vector3d(
			final(0,1), final(1,1), final(2,1));
		Eigen::Vector3d af = Eigen::Vector3d(
			final(0,2), final(1,2), final(2,2));

		Eigen::VectorXd delta(9);
		delta(0) = pf(0) - p0(0) - v0(0) * T - 0.5 * a0(0) * pow(T,2);
		delta(1) = pf(1) - p0(1) - v0(1) * T - 0.5 * a0(1) * pow(T,2);
		delta(2) = pf(2) - p0(2) - v0(2) * T - 0.5 * a0(2) * pow(T,2);
		delta(3) = (vf(0) - v0(0) - a0(0) * T);
		delta(4) = (vf(1) - v0(1) - a0(1) * T);
		delta(5) = (vf(2) - v0(2) - a0(2) * T);
		delta(6) = (af(0) - a0(0));
		delta(7) = (af(1) - a0(1));
		delta(8) = (af(2) - a0(2));

		Eigen::Matrix3d m;
		m << 720, -360*T, 60*pow(T,2),
			-360*T, 168*pow(T,2), -24*pow(T,3),
			60*pow(T,2), -24*pow(T,3), 3*pow(T,4);

        // Print out the delta vector
        // printf("delta:\n");
        // for (int i = 0; i < 9; i++)
        // {
        //     printf("%lf\n", delta(i));
        // }
        // printf("\n");

		Eigen::MatrixXd M(9,9);
		M.setZero();

		// Make the 3x3 matrix into a 9x9 matrix so that xyz is included
		for (int i = 0; i < 9; i++)
		{
			int l_m = (i+1) % 3 + 3*(int)(((i+1) % 3) == 0) - 1;
			int h_m = (int)ceil((double)(i+1) / 3.0) - 1;

			int l_M = (l_m)*3;
			int h_M = (h_m)*3;

			M(l_M,h_M) = m(l_m,h_m);
			M(l_M+1,h_M+1) = m(l_m,h_m);
			M(l_M+2,h_M+2) = m(l_m,h_m);

            // printf("M(%d,%d) to M(%d,%d)\n", l_M, l_M, l_M+2, h_M+2);
            // printf("m(%d,%d), %lf\n", l_m, h_m, m(l_m,h_m));
		}

        // Print out the M matrix vector
        // printf("M matrix:\n");
        // for (int i = 0; i < 9; i++)
        // {
        //     for (int j = 0; j < 9; j++)
        //     {
        //         printf("%lf ", M(i,j));
        //     }
        //     printf("\n");
        // }
        // printf("\n");

		Eigen::VectorXd abg(9);

		abg = (1/pow(T,5) * M * delta);
        *alpha = Eigen::Vector3d(abg(0), abg(1), abg(2));
        *beta = Eigen::Vector3d(abg(3), abg(4), abg(5));
        *gamma = Eigen::Vector3d(abg(6), abg(7), abg(8));

	}

    px4_array_container get_discrete_points(matrix::SquareMatrix<double, 3> initial,
		matrix::SquareMatrix<double, 3> final,
        double total_time, double command_time, matrix::Vector3d alpha, matrix::Vector3d beta,
		matrix::Vector3d gamma, int& waypoint_size)
	{
        // Update the initial values
		matrix::Vector3d p0 = matrix::Vector3d(
			initial(0,0), initial(1,0), initial(2,0));
		matrix::Vector3d v0 = matrix::Vector3d(
			initial(0,1), initial(1,1), initial(2,1));
		matrix::Vector3d a0 = matrix::Vector3d(
			initial(0,2), initial(1,2), initial(2,2));

		// Update the destination values
		matrix::Vector3d pf = matrix::Vector3d(
			final(0,0), final(1,0), final(2,0));
		matrix::Vector3d vf = matrix::Vector3d(
			final(0,1), final(1,1), final(2,1));
		matrix::Vector3d af = matrix::Vector3d(
			final(0,2), final(1,2), final(2,2));

        waypoint_size = (int)(total_time / command_time);
        // double corrected_interval = total_time / (double)waypoint_size;
        px4_array_container desired_states;
        for (int i = 0; i < waypoint_size; i++)
        {
            matrix::Vector3d pos = (
                alpha/120 * pow((command_time*i),5) + 
                beta/24 * pow((command_time*i),4) + 
                gamma/6 * pow((command_time*i),3) + 
                a0/2 * pow((command_time*i),2) + 
                v0 * command_time*i + 
				p0);
            matrix::Vector3d vel = (
				alpha/24 * pow((command_time*i),4) + 
                beta/6 * pow((command_time*i),3) + 
                gamma/2 * pow((command_time*i),2) + 
				a0 * command_time*i + 
				v0);
            matrix::Vector3d acc = (
				alpha/6 * pow((command_time*i),3) + 
                beta/2 * pow((command_time*i),2) + 
				gamma * command_time*i + 
				a0);
            
            double data[] = {
                pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2)
            };
            matrix::Vector<double, 9> desired_state(data);

            desired_states.push_back(desired_state);
        }
        return desired_states;
	}

	// check_z_vel using PX4 matrix (without eigen)
    int check_z_vel(matrix::SquareMatrix<double, 3> initial,
		matrix::SquareMatrix<double, 3> final,
        double total_time, double command_time, matrix::Vector3d alpha, matrix::Vector3d beta,
		matrix::Vector3d gamma)
	{
		matrix::Vector3d v0 = matrix::Vector3d(
			initial(0,1), initial(1,1), initial(2,1));
		matrix::Vector3d a0 = matrix::Vector3d(
			initial(0,2), initial(1,2), initial(2,2));

        int waypoint_size = (int)(total_time / command_time);
        double corrected_interval = total_time / (double)waypoint_size;
		int bad_counts = 0;
        for (int i = 0; i < waypoint_size; i++)
        {
            matrix::Vector3d vel = (alpha/24 * pow((corrected_interval*i),4) + 
                beta/6 * pow((corrected_interval*i),3) + 
                gamma/2 * pow((corrected_interval*i),2) + a0 * corrected_interval + v0);
            
            if (vel(2) > 0.001)
                bad_counts += 1;
        }

		return bad_counts;
	}

	// check_z_vel using eigen
	int check_z_vel(Eigen::Matrix3d initial,
		Eigen::Matrix3d final,
        double total_time, double command_time, Eigen::Vector3d alpha, Eigen::Vector3d beta,
		Eigen::Vector3d gamma)
	{
		Eigen::Vector3d v0 = Eigen::Vector3d(
			initial(0,1), initial(1,1), initial(2,1));
		Eigen::Vector3d a0 = Eigen::Vector3d(
			initial(0,2), initial(1,2), initial(2,2));

        int waypoint_size = (int)ceil(total_time / command_time);
        double corrected_interval = total_time / (double)waypoint_size;
		int bad_counts = 0;
        for (int i = 0; i < waypoint_size; i++)
        {
            Eigen::Vector3d vel = (alpha/24 * pow((corrected_interval*i),4) + 
                beta/6 * pow((corrected_interval*i),3) + 
                gamma/2 * pow((corrected_interval*i),2) + a0 * corrected_interval + v0);
            
            if (vel(2) > 0.001)
                bad_counts += 1;
        }

		return bad_counts;
	}
}

#endif