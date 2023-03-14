#pragma once

#include <string>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <vector>
#include <set>

#include "../header files/triangulation.h"

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Dynamic;

typedef typename Eigen::Triplet<double> T;

/*
* Calculates a reduced matrix from a given matrix by filtering by indices
* 
* @param M the full matrix
* @param fNodes a vector of integer indices
* @param dof the size of the fNodes vector
* @returns a submatrix of M that has the row and coloumn indices specified by fNodes
* 
* Example:
*	MatrixXd M {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
*	VectorXd fNodes {{0, 2}};
*	MatrixXd M_reduced = extractSparseMatrix(M, fNodes, 2);
*	cout << M_reduced << endl;
*	
*	1	3
*	3	9
*/
SparseMatrix<double> extractSparseMatrix(SparseMatrix<double> M, std::vector<int> fNodes, int dof);

/*
* Calculates a reduced matrix from a given matrix by filtering by indices
*
* @param M the full matrix
* @param xNodes a vector of column integer indices
* @param yNodes a vector of row integer indices
* @returns a submatrix of M that has the row indices specified by xNodes and column indices specified by yNodes
* 
* Example:
*	MatrixXd M {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
*	VectorXd xNodes {{0, 1}};
*	VectorXd yNodes {{1, 2}};
*	MatrixXd M_reduced = extractSparseMatrix(M, xNodes, yNodes);
*	cout << M_reduced << endl;
*
*	4	5
*	7	8
*/
SparseMatrix<double> extractSparseMatrix(SparseMatrix<double> M, std::vector<int> xNodes, std::vector<int> yNodes);

/*
* Calculates a reduced matrix from a given matrix by filtering by column indices
*
* @param M the full matrix
* @param xNodes a vector of column integer indices
* @param yNodes a vector of row integer indices
* @returns a submatrix of M that has the column indices specified by yNodes
*
* Example:
*	MatrixXd M {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
*	VectorXd yNodes {{1, 2}};
*	MatrixXd M_reduced = extractSparseMatrix(M, xNodes, yNodes);
*	cout << M_reduced << endl;
*
*	4	5	6
*	9	8	9
*/
SparseMatrix<double> extractVerticalSparseMatrix(SparseMatrix<double> M, std::vector<int> nodes, int dof);

/*
* Calculates a vector that is filled with a specified value for missing indices in fNodes
*
* @param V the specified vector with the same size as fNodes
* @param nodes indices that are present in V from the perspective of "V_full"
* @param full the size of V_full
* @param value the value to be filled
* @returns returns V_full that has the corresponding value of V if index is in nodes else value
*
* Example:
*	VectorXd V = {{10, 20, 30, 40, 50, 60, 80}}
*	VectorXd nodes {{1, 2, 4, 7, 11, 12, 13}};
*	int full = 14
*	int value = -1
*	VectorXd V_full = insertDenseVector(V, fNodes, full, value);
*	cout << V_full << endl;
*
*	-1
*	10
*	20
*	-1
*	30
*	-1
*	-1
*	40
*	-1
*	-1
*	-1
*	50
*	60
*	80
*	-1
*/
VectorXd insertDenseVector(VectorXd V, std::vector<int> fNodes, int full, double value);

class P1FiniteElement {
public:
	// Stores the triangulation data
	Triangulation tri;

	// stores the length/areas/volumes of the simplicies defined by tri
	// has size tri.num_edges
	Vector<double, Dynamic> volumes;
	
	// stores the length/areas of the simplicies defined by tri on the Neumann boundary
	// has size tri.Nb.size()
	Vector<double, Dynamic> volumes_Nb;
	
	// stores the barycenters of the simplices defined by tri
	// the matrix has dimensions tri.node_dim x tri.num_edges
	// the i-th column stores the barycenter of the i-th simplice in the triangulation
	MatrixXd bctr;

	// stores the barycenters of the edges defined by tri on the Neumann bounday
	// the matrix has dimensions tri.node_dim x tri.Nb.size()
	// the i-th column stores the barycenter of the i-th edge in the triangulation on the Neumann boundary
	MatrixXd bctr_Nb;

	// degrees of freedom
	int dof = 0;

	// ordered vector of free nodes
	std::vector<int> fNodes;

	// set of free nodes
	std::vector<int> ifNodes;

	// array of matricies that store the gradients of the nodal basis functions of the triangulation
	// the length of the array is tri.num_nodes. Each matrix has dimension tri.node_dim x tri.edge_dim
	MatrixXd* grads;

	SparseMatrix<double> s;
	SparseMatrix<double> sff;
	SparseMatrix<double> sr;

	SparseMatrix<double> m;
	SparseMatrix<double> mff;
	SparseMatrix<double> mr;

	Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower || Eigen::Upper> solver;

public:
	// sets up a finite element given a triangulation
	P1FiniteElement(Triangulation tri);

	//
	//void initStiffness();

	// returns the stiffness matrix for the finite element
	SparseMatrix<double> getStiffness();

	// returns the mass matrix for the finite element
	SparseMatrix<double> getMass();

	//VectorXd getDiffusionIntegrator(std::function<MatrixXd> A);

	/*
	* Returns a Vector for the FEM linear system corresponding to a right hand side 
	* mass integration of fixed function. The vector is created by zero-order integration
	* i.e. on each simplice, f is approximated by f(m_s) where m_s is the barycenter of
	* said simplice
	* 
	* @param f a analytical function from which the Vector is calculated
	* 
	* Explanation:
	*	regard -laplace(u) = f
	*	then its weak formulation is (assuming Dirichlet boundary):
	*		integral grad(u) * grad(v) = integral f * v
	*	the return value of this function is the FEM Vector for the right hand side
	*/
	VectorXd getRHSMassIntegrator(std::function<double(VectorXd&)> f);

	/*
	* Returns a Vector for the FEM linear system corresponding to the Neumann boundary given 
	* a fixed function g. The vector is created by zero-order integration
	* i.e. on each Neumann boundary edge, g is approximated by g(m_s) where m_s is the barycenter of
	* said Neumann boundary edge
	*
	* @param g a analytical function from which the Vector is calculated
	*
	* Explanation:
	*	regard -laplace(u) = f with neumann boundary g
	*	then its weak formulation is
	*		integral grad(u) * grad(v) = integral f * v + integral g * v
	*	where the latter integral is over the boundary of the domain
	*	the return value of this function is the FEM Vector last integral
	*/
	VectorXd getNeumannBoundaryMassIntegrator(std::function<double(VectorXd&)> g);

	/*
	* creates a c4n gridfunction from a mathematical function
	*/
	VectorXd nodal_interpolant(std::function<double(VectorXd)> f);

private:
	void getStiffnessThread(SparseMatrix<double>& s, SparseMatrix<double>& sff,
							SparseMatrix<double>& sr, int index, int num_threads);

	void getMassThread(SparseMatrix<double>& s, SparseMatrix<double>& sff,
		SparseMatrix<double>& sr, int index, int num_threads);

	void InitP1FiniteElementThread(Triangulation& t, int index, int num_threads);
};