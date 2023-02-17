#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <map>
#include <set>

#include "../header files/triangulation.h"

using Eigen::Vector;
using Eigen::Matrix;
using Eigen::SparseMatrix;

Triangulation::Triangulation(std::string path) {
	std::string line;
	int num_nodes = 0;
	int num_edges = 0;
	int node_dim = 0;
	int edge_dim = 0;
	Matrix<double, Dynamic, Dynamic> c4n;
	Matrix<int, Dynamic, Dynamic> n4e;
	Matrix<int, Dynamic, Dynamic> Db;

	std::ifstream infile("test_file.txt");

	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		std::string key;
		iss >> key;

		//std::cout << key << std::endl;

		if (key == "node_dim") {
			iss >> node_dim;
		}
		else if (key == "num_nodes") {
			iss >> num_nodes;
			c4n.resize(num_nodes, node_dim);
		}
		else if (key == "c4n") {
			for (int i = 0; i < num_nodes; i++) {
				for (int j = 0; j < node_dim; j++) {
					iss >> c4n(i, j);
				}
			}
		}
		else if (key == "edge_dim") {
			iss >> edge_dim;
		}
		else if (key == "num_edges") {
			iss >> num_edges;
			n4e.resize(num_edges, edge_dim);
		}
		else if (key == "n4e") {
			for (int i = 0; i < num_edges; i++) {
				for (int j = 0; j < 3; j++) {
					iss >> n4e(i, j);
				}
			}
		}
		else if (key == "db") {
			Db.resize(num_nodes, node_dim);
			for (int i = 0; i < num_nodes; i++) {
				for (int j = 0; j < 2; j++) {
					iss >> Db(i, j);
				}
			}
		}
	}

	this->node_dim = node_dim;
	this->num_nodes = num_nodes;
	this->edge_dim = edge_dim;
	this->num_edges = num_edges;
	this->c4n = c4n;
	this->n4e = n4e;
	this->n4e = Db;
}


Triangulation::Triangulation() {
	this->num_nodes = 4;
	this->num_edges = 2;
	this->node_dim = 2;
	this->edge_dim = 3;

	Matrix<double, 4, 2> c4n{ {0, 0}, {1, 0}, {1, 1}, {0, 1} };
	Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
	Matrix<int, 4, 2> Db{ {0, 1}, {1, 2}, {2, 3}, {3, 0} };
	//Matrix<int, 0, 0> Nb{ {} };

	this->c4n = c4n;
	this->n4e = n4e;
	this->Db = Db;
	//this->Nb = Nb;
}

void Triangulation::uniform_refinement() {
	std::map<std::pair<int,int>, int> mpts;
	std::vector<VectorXd>c4n;
	std::vector<Vector<int, Dynamic>> n4e;
	std::vector<Vector<int, Dynamic>> Db;

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

	for (int i = 0; i < this->Db.rows(); i++) {
		Vector<int, Dynamic> v1{ {this->Db(i, 0), mpts[{this->Db(i, 0), this->Db(i, 1)}]} };
		Vector<int, Dynamic> v2{ {mpts[{this->Db(i, 0), this->Db(i, 1)}], this->Db(i, 1)} };
		Db.push_back(v1);
		Db.push_back(v2);
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

	this->num_edges = num_edges_red;
	this->num_nodes = num_nodes_red;
	this->c4n = c4n_red;
	this->n4e = n4e_red;
	this->Db = Db_red;
}

void Triangulation::export_tri(std::string name) {
	// Open the output file
	std::ofstream outfile_c4n("meshes/" + name + "_c4n.tsv");
	if (!outfile_c4n.is_open()) {
		std::cerr << "Failed to open output file c4n." << std::endl;
		return;
	}

	// Write the matrix to the output file in TSV format
	for (int i = 0; i < c4n.rows(); i++) {
		for (int j = 0; j < c4n.cols(); j++) {
			outfile_c4n << c4n(i, j);
			if (j < c4n.cols() - 1) {
				outfile_c4n << "\t"; // use tab separator
			}
		}
		outfile_c4n << std::endl;
	}

	// Close the output file
	outfile_c4n.close();

	// Open the output file
	std::ofstream outfile_n43("meshes/" + name + "_n4e.tsv");
	if (!outfile_n43.is_open()) {
		std::cerr << "Failed to open output file n4e." << std::endl;
		return;
	}

	// Write the matrix to the output file in TSV format
	for (int i = 0; i < n4e.rows(); i++) {
		for (int j = 0; j < n4e.cols(); j++) {
			outfile_n43 << n4e(i, j);
			if (j < n4e.cols() - 1) {
				outfile_n43 << "\t"; // use tab separator
			}
		}
		outfile_n43 << std::endl;
	}

	// Close the output file
	outfile_n43.close();

	// Open the output file
	std::ofstream outfile_db("meshes/" + name + "_db.tsv");
	if (!outfile_db.is_open()) {
		std::cerr << "Failed to open output file db." << std::endl;
		return;
	}

	// Write the matrix to the output file in TSV format
	for (int i = 0; i < Db.rows(); i++) {
		for (int j = 0; j < Db.cols(); j++) {
			outfile_db << Db(i, j);
			if (j < Db.cols() - 1) {
				outfile_db << "\t"; // use tab separator
			}
		}
		outfile_db << std::endl;
	}

	// Close the output file
	outfile_db.close();

	std::cout << "Successfully exported " + name << std::endl;
}