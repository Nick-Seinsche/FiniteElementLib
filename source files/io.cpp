#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include "../header files/io.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;

std::ofstream prepareFile(IO::ExportType et, std::string name) {
	std::string path = "";
	if (et == IO::SOLUTION) path = "solutions/";
	if (et == IO::TRIANGULATION) path = "meshes/";

	std::ofstream outfile(path + name + ".tsv");
	if (!outfile.is_open()) {
		std::cerr << "Failed to export " + name << std::endl;
	}

	return outfile;
}

template <typename Matrix>
void IO::export_matrix_as_tsv(Matrix M, ExportType et, std::string name) {
	if (et == IO::DISABLED) return;
	std::ofstream outfile = prepareFile(et, name);

	// Write the matrix to the output file in TSV format
	for (int i = 0; i < M.rows(); i++) {
		for (int j = 0; j < M.cols(); j++) {
			outfile << M(i, j);
			if (j < M.cols() - 1) {
				outfile << "\t"; // use tab separator
			}
		}
		outfile << std::endl;
	}
}

template void IO::export_matrix_as_tsv<MatrixXd>(MatrixXd M, ExportType et, std::string name);
template void IO::export_matrix_as_tsv<MatrixXi>(MatrixXi M, ExportType et, std::string name);


template <typename Vector>
void IO::export_vector_as_tsv(Vector V, ExportType et, std::string name) {
	if (et == IO::DISABLED) return;
	std::ofstream outfile = prepareFile(et, name);

	// Write the matrix to the output file in TSV format
	for (int i = 0; i < V.rows(); i++) {
		outfile << V(i) << std::endl;
	}
}

template void IO::export_vector_as_tsv<VectorXd>(VectorXd V, ExportType et, std::string name);


MatrixXd import_tsv_as_matrix(std::string path, std::string name) {
	std::ifstream infile(path + name + ".tsv");
	std::vector<std::vector<double>> data;
	while (!infile.eof()) {
		std::string str;
		std::getline(infile, str);
		std::stringstream buffer(str);
		std::string temp;
		std::vector<double> values;

		while (getline(buffer, temp, '\t')) {
			values.push_back(::strtod(temp.c_str(), 0));
		}

		data.push_back(values);
	}

	int m = data.size();
	int n = data[0].size();
	MatrixXd out(m, n);

	for (int i = 0; i < data.size(); ++i) {
		for (int j = 0; j < data[0].size(); ++j) {
			out(i, j) = data[i][j];
		}
	}

	return out;
}