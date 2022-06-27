#include <iostream>
#include <string>
#include "obvp.h"

using namespace obvp;

/**
 * @brief bvp test case
 */

int main(int argc, char **argv) 
{
    matrix::Vector3d alpha, beta, gamma;
    double total_time = 5.0;
    matrix::SquareMatrix<double, 3> initial_state;
    initial_state.setZero();
    // Position
    initial_state.col(0) = matrix::Vector3d{5, -5, 10};
    // Velocity
    initial_state.col(1) = matrix::Vector3d{1, -1, -2};
    // Acceleration
    initial_state.col(2) = matrix::Vector3d{0, 0, 0};

    matrix::SquareMatrix<double, 3> final_state;
    final_state.setZero();
    // Position
    final_state.col(0) = matrix::Vector3d{0, 0, 0.5};
    // Velocity
    final_state.col(1) = matrix::Vector3d{0, 0, 0};
    // Acceleration
    final_state.col(2) = matrix::Vector3d{0, 0, 0};

    get_bvp_coefficients(initial_state, final_state, total_time,
		&alpha, &beta, &gamma);

    printf("alpha coefficient(%lf, %lf, %lf)\n", alpha(0), alpha(1), alpha(2));
    printf("beta coefficient(%lf, %lf, %lf)\n", beta(0), beta(1), beta(2));
    printf("gamma coefficient(%lf, %lf, %lf)\n", gamma(0), gamma(1), gamma(2));

    return 0;
}