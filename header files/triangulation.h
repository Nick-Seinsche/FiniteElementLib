#pragma once

#include <string>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Dynamic;

/*
* Storage class for Triangulations and Boundaries that also handles
* uniform refinement and i/o.
*/
class Triangulation {

public:
	// dimension of the nodes
	int node_dim;

	// dimension of the edges
	int edge_dim;
	
	// number of nodes
	int num_nodes;

	// number of deges
	int num_edges;

	// coordinates for nodes
	MatrixXd c4n;

	// nodes for edges
	Matrix<int, Dynamic, Dynamic> n4e;

	// Dirichlet boundary
	Matrix<int,Dynamic, Dynamic> Db;

	// Neumann boundary
	Matrix<int, Dynamic, Dynamic> Nb;

	/*
	* Default triangulation
	*/
	Triangulation();

	/*
	* Pass a triangulation as argument
	*/
	Triangulation(MatrixXd c4n, MatrixXi n4e, MatrixXi Db, MatrixXi Nb);

	/*
	* Loads triangulation data from a file
	* @param path the path of the file
	*/
	Triangulation(std::string path);

	/*
	* Updates the triangulation with the uniform refined 
	* triangulation data
	*/
	void uniform_refinement();

	/*
	* Exports the triangulation data to the meshes directory
	* 
	* @param the name of the triangulation data files that are
	*	created e.g. <name>_c4n.tsv
	*/
	void export_tri(std::string name);
};

