#pragma once

#include <Eigen/dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

void export_func(VectorXd sol, std::string name);

VectorXd create(std::function<double(VectorXd)> f, MatrixXd grid, int num_nodes);