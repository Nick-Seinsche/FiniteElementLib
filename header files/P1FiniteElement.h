#pragma once

#include <string>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
#include <set>

#include "../header files/triangulation.h"

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Dynamic;

typedef Eigen::Triplet<double> T;

SparseMatrix<double> extractSparseMatrix(SparseMatrix<double> M, std::vector<int> fNodes, int dof);
SparseMatrix<double> extractSparseMatrix(SparseMatrix<double> M, std::vector<int> xNodes, std::vector<int> yNodes);
SparseMatrix<double> extractVerticalSparseMatrix(SparseMatrix<double> M, std::vector<int> nodes);
VectorXd extractDenseVector(VectorXd V, std::vector<int> fNodes, int dof);
VectorXd insertDenseVector(VectorXd V, std::vector<int> fNodes, int full);

class P1FiniteElement {
public:
	Triangulation tri;

	Vector<double, Dynamic> volumes;
	MatrixXd bctr;

	int dof = 0;
	std::vector<int> fNodes;
	std::set<int> sfNodes;

	MatrixXd* grads;
	SparseMatrix<double> s;
	SparseMatrix<double> m;
	VectorXd rhs;

public:
	P1FiniteElement(Triangulation tri);

	SparseMatrix<double> getStiffness();

	SparseMatrix<double> getFreeNodeStiffness();

	SparseMatrix<double> getMass();

	SparseMatrix<double> getFreeNodeMass();

	//VectorXd getDiffusionIntegrator(std::function<MatrixXd> A);

	VectorXd getRHSMassIntegrator(std::function<double(VectorXd&)> f);

	VectorXd getFreeNodeRHSMassIntegrator(std::function<double(VectorXd&)> f);

	VectorXd addDirichletNodes(VectorXd X);
};