#pragma once

#include <vector>
#include <Eigen/Dense>
#include <iostream>
//#include "../header files/triangulation.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef struct {
	std::vector<MatrixXd> points;
	std::vector<VectorXd> weights;

	// todo
	std::vector<MatrixXd> pointsBoundary;
	std::vector<VectorXd> weightsBoundary;
} Integrator;


namespace integrators {
	//Integrator zeroOrderMidpoint(Triangulation tri);

	//Integrator gauss3Point(Triangulation tri);
}

