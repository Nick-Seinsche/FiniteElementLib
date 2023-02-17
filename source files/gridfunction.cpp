#include <iostream>
#include <fstream>
#include "../header files/gridfunction.h"

void export_func(VectorXd sol, std::string name) {

	// Open the output file
	std::ofstream outfile("solutions/" + name + ".tsv");
	if (!outfile.is_open()) {
		std::cerr << "Failed to open output file." << std::endl;
		return;
	}

	// Write the matrix to the output file in TSV format
	for (int i = 0; i < sol.size(); i++) {
			outfile << sol(i);
		outfile << std::endl;
	}

	// Close the output file
	outfile.close();

	std::cout << "Successfully exported " + name << std::endl;
}

VectorXd create(std::function<double(VectorXd)> f, MatrixXd grid, int num_nodes) {
	VectorXd gridf(num_nodes);
	for (int i = 0; i < num_nodes; i++) {
		gridf(i) = f(grid.row(i).transpose());
	}
	return gridf;
}