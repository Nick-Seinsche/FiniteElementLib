#pragma once

#include <string>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
#include <set>

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Dynamic;

class Triangulation {

public:
	int node_dim;
	int edge_dim;
	int num_nodes;
	int num_edges;

	MatrixXd c4n;
	Matrix<int, Dynamic, Dynamic> n4e;
	Matrix<int,Dynamic, Dynamic> Db;
	Matrix<int, Dynamic, Dynamic> Nb;

	Triangulation();

	Triangulation(std::string path);

	void uniform_refinement();

	void export_tri(std::string name);
};

