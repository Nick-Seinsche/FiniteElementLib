#include "../header files/integrators.h"

/*
Integrator integrators::gauss3Point(Triangulation tri) {
	Integrator ig;

	MatrixXd Q_T(tri.node_dim, tri.edge_dim - 1);
	VectorXd gaussweight(tri.edge_dim);
	VectorXd gausspoint(tri.edge_dim);

	VectorXd unit_gaussweights(3);
	MatrixXd unit_gausspoints;

	if (tri.node_dim == 2 && tri.edge_dim == 3) {
		unit_gausspoints(tri.node_dim, 3);
		unit_gausspoints << 1 / 6, 4 / 6, 1 / 6,
							1 / 6, 1 / 6, 4 / 6;
		unit_gaussweights << 1 / 6, 1 / 6, 1 / 6;
	}
	else if (tri.node_dim == 1 && tri.edge_dim == 2) {
		std::cerr << "not implemented" << std::endl;
	}

	for (int i = 0; i < tri.num_edges; i++) {
		for (int j = 0; j < tri.edge_dim - 1; j++) {
			Q_T.col(j) = (tri.c4n.row(tri.n4e(i, j + 1)) - tri.c4n.row(tri.n4e(i, 0))).transpose();
		}
		for (int j = 0; j < tri.edge_dim; j++) {
			gausspoint = tri.c4n.row(tri.n4e(i, 0)).transpose() + Q_T * unit_gausspoints.col(j);
			gaussweight(j) = std::abs(Q_T.determinant()) * unit_gaussweights(j);
		}
		ig.weights.push_back(gaussweight);
		ig.points.push_back(gausspoint);
	}

	return ig;
}

Integrator integrators::zeroOrderMidpoint(Triangulation tri) {
	Integrator ig;
	MatrixXd bctr(tri.node_dim, tri.num_edges);
	
	for (int i = 0; i < tri.num_edges; i++) {
		for (int j = 0; j < tri.edge_dim; j++) {
			bctri += tri.c4n.row(tri.n4e(i, j)).transpose();
		}
		bctri = VectorXd::Zero(tri.node_dim);
		bctr.col(i) = bctri / tri.edge_dim;
	}
	
	VectorXd bctrnbi;

	for (int i = 0; i < tri.Nb.rows(); i++) {
		bctrnbi = VectorXd::Zero(tri.node_dim);
		for (int j = 0; j < tri.edge_dim - 1; j++) {
			bctrnbi += tri.c4n.row(tri.Nb(i, j)).transpose();
		}
		bctr_Nb.col(i) = bctrnbi / (tri.edge_dim - 1);
	}
	MatrixXd bctr_Nb(tri.node_dim, tri.Nb.rows());
	VectorXd bctri;

	
}
*/