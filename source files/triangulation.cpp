#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <map>

#include "../header files/triangulation.h"
#include "../header files/io.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SparseMatrix;

Triangulation::Triangulation(std::string name) {
	std::ifstream file_c4n("meshes/" + name + "_c4n.tsv");
	std::string line;
	std::vector<std::vector<double>> data;

	// Read the file line by line
	while (std::getline(file_c4n, line)) {
		std::vector<double> row_data;
		std::string field;

		// Parse each line by tab delimiter
		std::istringstream line_stream(line);
		while (std::getline(line_stream, field, '\t')) {
			row_data.push_back(std::stod(field));
		}
		data.push_back(row_data);
	}

	// Create an Eigen MatrixXd from the parsed data
	int num_nodes = static_cast<int>(data.size());
	int node_dim = static_cast<int>(data[0].size());
	MatrixXd c4n(num_nodes, node_dim);
	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < node_dim; j++) {
			c4n(i, j) = data[i][j];
		}
	}

	data.empty();

	std::ifstream file_n4e("meshes/" + name + "_n4e.tsv");

	// Read the file line by line
	while (std::getline(file_n4e, line)) {
		std::vector<double> row_data;
		std::string field;

		// Parse each line by tab delimiter
		std::istringstream line_stream(line);
		while (std::getline(line_stream, field, '\t')) {
			row_data.push_back(std::stod(field));
		}
		data.push_back(row_data);
	}

	// Create an Eigen MatrixXd from the parsed data
	int num_edges = static_cast<int>(data.size());
	int edge_dim = static_cast<int>(data[0].size());
	Matrix<int, Dynamic, Dynamic> n4e(num_edges, edge_dim);
	for (int i = 0; i < num_edges; i++) {
		for (int j = 0; j < edge_dim; j++) {
			n4e(i, j) = data[i][j];
		}
	}

	this->node_dim = node_dim;
	this->num_nodes = num_nodes;
	this->c4n = c4n;
	this->edge_dim = edge_dim;
	this->num_edges = num_edges;
	this->n4e = n4e;
	
	std::cerr << "todo" << std::endl;

	//this->Db = Db;
	//this->Nb = Nb;
}


Triangulation::Triangulation() {
	this->num_nodes = 4;
	this->num_edges = 2;
	this->node_dim = 2;
	this->edge_dim = 3;

	Matrix<double, 4, 2> c4n{ {0, 0}, {1, 0}, {1, 1}, {0, 1} };
	Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
	Matrix<int, 1, 2> Db{ {1, 2} };
	Matrix<int, 3, 2> Nb{ {0, 1}, {2, 3}, {3, 0} };

	this->c4n = c4n;
	this->n4e = n4e;
	this->Db = Db;
	this->Nb = Nb;
}

Triangulation::Triangulation(MatrixXd c4n, MatrixXi n4e, MatrixXi Db, MatrixXi Nb) {
	this->num_nodes = static_cast<int>(c4n.rows());
	this->num_edges = static_cast<int>(n4e.rows());
	this->node_dim = static_cast<int>(c4n.cols());
	this->edge_dim = static_cast<int>(n4e.cols());

	this->c4n = c4n;
	this->n4e = n4e;
	this->Db = Db;
	this->Nb = Nb;
}

void Triangulation::uniform_refinement() {
	std::map<std::pair<int,int>, int> mpts;
	std::vector<VectorXd>c4n;
	std::vector<Vector<int, Dynamic>> n4e;
	std::vector<Vector<int, Dynamic>> Db;
	std::vector<Vector<int, Dynamic>> Nb;

	for (int i = 0; i < num_nodes; i++) {
		c4n.push_back(this->c4n.row(i));
	}

	int num_nodes_red = num_nodes;
	int num_edges_red = 0;
	for (int i = 0; i < num_edges; i++) {
		for (int j = 0; j < edge_dim; j++) {
			for (int k = j + 1; k < edge_dim; k++) {
				//n4e(i, j) n4e(i, k);
				if (mpts.find({ this->n4e(i, j), this->n4e(i, k)}) == mpts.end()) {
					c4n.push_back(0.5 * (c4n[this->n4e(i, j)] + c4n[this->n4e(i, k)]));
					mpts[{ this->n4e(i, j), this->n4e(i, k) }] = num_nodes_red;
					mpts[{ this->n4e(i, k), this->n4e(i, j) }] = num_nodes_red;
					num_nodes_red++;
				}
			}
		}
	}

	for (int i = 0; i < num_edges; i++) {
		if (node_dim == 1) {
			Vector<int, Dynamic> v1 { {this->n4e(i, 0), mpts[{ this->n4e(i, 0), this->n4e(i, 1) }]} };
			Vector<int, Dynamic> v2 { {mpts[{ this->n4e(i, 0), this->n4e(i, 1)}], this->n4e(i, 1)} };
			
			n4e.push_back(v1);
			n4e.push_back(v2);
			
			num_edges_red += 2;
		}
		else if (node_dim == 2) {
			Vector<int, Dynamic> v1{ {this->n4e(i, 0), mpts[{ this->n4e(i, 0), this->n4e(i, 1) }], mpts[{ this->n4e(i, 0), this->n4e(i, 2) }]}};
			Vector<int, Dynamic> v2{ {this->n4e(i, 1), mpts[{ this->n4e(i, 1), this->n4e(i, 2) }], mpts[{ this->n4e(i, 1), this->n4e(i, 0) }]} };
			Vector<int, Dynamic> v3{ {this->n4e(i, 2), mpts[{ this->n4e(i, 2), this->n4e(i, 0) }], mpts[{ this->n4e(i, 2), this->n4e(i, 1) }]} };
			Vector<int, Dynamic> v4{ {mpts[{ this->n4e(i, 0),  this->n4e(i, 1) }], mpts[{  this->n4e(i, 1), this->n4e(i, 2) }], mpts[{ this->n4e(i, 2), this->n4e(i, 0) }]} };

			n4e.push_back(v1);
			n4e.push_back(v2);
			n4e.push_back(v3);
			n4e.push_back(v4);

			num_edges_red += 4;
		}
	}

	if (node_dim == 1) {
		std::cerr << "not implemented" << std::endl;
	} else if (node_dim == 2) {
		for (int i = 0; i < this->Db.rows(); i++) {
			Vector<int, Dynamic> v1{ {this->Db(i, 0), mpts[{this->Db(i, 0), this->Db(i, 1)}]} };
			Vector<int, Dynamic> v2{ {mpts[{this->Db(i, 0), this->Db(i, 1)}], this->Db(i, 1)} };
			Db.push_back(v1);
			Db.push_back(v2);
		}

		for (int i = 0; i < this->Nb.rows(); i++) {
			Vector<int, Dynamic> v1{ {this->Nb(i, 0), mpts[{this->Nb(i, 0), this->Nb(i, 1)}]} };
			Vector<int, Dynamic> v2{ {mpts[{this->Nb(i, 0), this->Nb(i, 1)}], this->Nb(i, 1)} };
			Nb.push_back(v1);
			Nb.push_back(v2);
		}
	}
	else {
		std::cerr << "not implemented" << std::endl;
	}

	MatrixXd c4n_red(num_nodes_red, node_dim);

	for (int i = 0; i < num_nodes_red; i++) {
		c4n_red.row(i) = c4n[i].transpose();
	}

	Matrix<int, Dynamic, Dynamic> n4e_red(num_edges_red, edge_dim);

	for (int i = 0; i < num_edges_red; i++) {
		n4e_red.row(i) = n4e[i].transpose();
	}

	Matrix<int, Dynamic, Dynamic> Db_red(Db.size(), edge_dim - 1);

	for (int i = 0; i < Db.size(); i++) {
		Db_red.row(i) = Db[i].transpose();
	}

	Matrix<int, Dynamic, Dynamic> Nb_red(Nb.size(), edge_dim - 1);

	for (int i = 0; i < Nb.size(); i++) {
		Nb_red.row(i) = Nb[i].transpose();
	}

	this->num_edges = num_edges_red;
	this->num_nodes = num_nodes_red;
	this->c4n = c4n_red;
	this->n4e = n4e_red;
	this->Db = Db_red;
	this->Nb = Nb_red;
}

void Triangulation::export_tri(std::string name) {
	// Open the output file
	IO::export_matrix_as_tsv<MatrixXd>(c4n, IO::TRIANGULATION, name + "_c4n");
	
	IO::export_matrix_as_tsv<MatrixXi>(n4e, IO::TRIANGULATION, name + "_n4e");

	IO::export_matrix_as_tsv<MatrixXi>(Db, IO::TRIANGULATION, name + "_db");

	IO::export_matrix_as_tsv<MatrixXi>(Nb, IO::TRIANGULATION, name + "_nb");

	std::cout << "Successfully exported " + name << std::endl;
}