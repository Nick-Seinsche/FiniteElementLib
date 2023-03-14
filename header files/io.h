#pragma once

#include <Eigen/Dense>

namespace IO {
	enum ExportType {SOLUTION, TRIANGULATION, DISABLED};

	// writes Eigen::Matrix M to specified file in tsv format
	template <typename Matrix>
	void export_matrix_as_tsv(Matrix M, ExportType et, std::string name);

	// writes Eigen::Vector V to specified file in tsv format
	template <typename Vector>
	void export_vector_as_tsv(Vector V, ExportType et, std::string name);
}