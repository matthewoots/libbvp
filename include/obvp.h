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
#include "math.hpp"
#include "geo.h"
#include "mathlib.h"

using namespace std;
using matrix::Vector2d;
using matrix::Vector2f;
using matrix::wrap_pi;
using matrix::wrap_2pi;

namespace obvp
{
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
        //     printf("%lf\n", delta(i));
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
        *alpha = matrix::Vector3d(abg(0), abg(1), abg(2));
        *beta = matrix::Vector3d(abg(3), abg(4), abg(5));
        *gamma = matrix::Vector3d(abg(6), abg(7), abg(8));

	}
}

#endif