#include <iostream>

#include "../header files/P1FiniteElement.h"

int __factorial(int k) {
	int kf = 1;
	for (int i = 2; i <= k; ++i) {
		kf *= i;
	}
	return kf;
}

SparseMatrix<double> extractSparseMatrix(SparseMatrix<double> M, std::vector<int> fNodes, int dof) {
	SparseMatrix<double> M_red(dof, dof);
	for (int k = 0; k < M.outerSize(); ++k)
		for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
			auto rw = std::find(fNodes.begin(), fNodes.end(), it.col());
			auto rc = std::find(fNodes.begin(), fNodes.end(), it.row());
			if (rw != fNodes.end() && rc != fNodes.end()) {
				M_red.insert(rc - fNodes.begin(), rw - fNodes.begin()) = it.value();
			}
		}
	return M_red;
}

SparseMatrix<double> extractSparseMatrix(SparseMatrix<double> M, std::vector<int> xNodes, std::vector<int> yNodes) {
	int nx = xNodes.size();
	int ny = yNodes.size();
	SparseMatrix<double> M_red(ny, nx);
	for (int k = 0; k < M.outerSize(); ++k)
		for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
			auto rw = std::find(xNodes.begin(), xNodes.end(), it.col());
			auto rc = std::find(yNodes.begin(), yNodes.end(), it.row());
			if (rw != xNodes.end() && rc != yNodes.end()) {
				M_red.insert(rc - yNodes.begin(), rw - xNodes.begin()) = it.value();
			}
		}
	return M_red;
}

SparseMatrix<double> extractVerticalSparseMatrix(SparseMatrix<double> M, std::vector<int> nodes) {
	int nx = nodes.size();
	SparseMatrix<double> M_red(nx, M.cols());
	for (int k = 0; k < M.outerSize(); ++k)
		for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
			auto rw = std::find(nodes.begin(), nodes.end(), it.col());
			if (rw != nodes.end()) {
				M_red.insert(rw - nodes.begin(), it.col()) = it.value();
			}
		}
	return M_red;
}

VectorXd extractDenseVector(VectorXd V, std::vector<int> fNodes, int dof) {
	VectorXd V_red(dof);
	int j = 0;
	for (int i = 0; i < V.size(); i++) {
		if (std::find(fNodes.begin(), fNodes.end(), i) != fNodes.end()) {
			V_red[j] = V[i];
			j++;
		}
	}
	return V_red;
}

VectorXd insertDenseVector(VectorXd V, std::vector<int> fNodes, int full) {
	VectorXd V_full(full);
	int j = 0;
	for (int i = 0; i < full; i++) {
		if (std::find(fNodes.begin(), fNodes.end(), i) != fNodes.end()) {
			V_full(i) = V(j);
			j++;
		} else {
			V_full(i) = 0;
		}
	}
	return V_full;
}

P1FiniteElement::P1FiniteElement(Triangulation tri) {
	this->tri = tri;

	MatrixXd X_T(tri.node_dim + 1, tri.node_dim + 1);
	MatrixXd* grads = new MatrixXd[tri.num_edges];
	VectorXd vol(tri.num_edges);

	MatrixXd bctr(tri.node_dim, tri.num_edges);

	std::set<int> sdb;
	for (auto k : tri.Db.reshaped()) {
		sdb.insert(k);
	}
	std::vector<int> fNodes;
	std::set<int> sfNodes;
	for (int i = 0; i < tri.num_nodes; i++) {
		if (sdb.find(i) == sdb.end()) {
			fNodes.push_back(i);
			sfNodes.insert(i);
			this->dof++;
		}
	}

	VectorXd rhs = VectorXd::Zero(tri.num_nodes);

	MatrixXd I_d(tri.node_dim + 1, tri.node_dim);
	I_d.row(0) = MatrixXd::Zero(1, tri.node_dim);
	I_d.block(1, 0, tri.node_dim, tri.node_dim)
		= MatrixXd::Identity(tri.node_dim, tri.node_dim);

	VectorXd bctri;

	for (int i = 0; i < tri.num_edges; i++) {
		bctri = VectorXd::Zero(tri.node_dim);


		for (int j = 0; j < tri.edge_dim; j++) {
			X_T(0, j) = 1;
			X_T.block(1, j, tri.node_dim, 1) = tri.c4n.row(tri.n4e(i, j)).transpose();
			bctri += tri.c4n.row(tri.n4e(i, j)).transpose();
		}

		vol(i) = X_T.determinant() / __factorial(tri.node_dim);

		bctr.block(0, i, tri.node_dim, 1) = bctri / tri.edge_dim;

		grads[i] = (X_T.inverse() * I_d).transpose();
	}

	this->volumes = vol;
	this->grads = grads;
	this->bctr = bctr;
	this->rhs = rhs;
	this->fNodes = fNodes;
	this->sfNodes = sfNodes;
}

SparseMatrix<double> P1FiniteElement::getStiffness() {
	std::vector<T> tripletList;
	double value = 0;

	for (int t = 0; t < tri.num_edges; t++) {
		for (int i = 0; i < tri.edge_dim; i++) {
			for (int j = i; j < tri.edge_dim; j++) {
				value = volumes(t) * (grads[t].col(i).transpose() * grads[t].col(j)).value();

				tripletList.push_back(T(tri.n4e(t, i), tri.n4e(t, j), value));
				if (i != j) tripletList.push_back(T(tri.n4e(t, j), tri.n4e(t, i), value));
			}
		}
	}
	SparseMatrix<double> s(tri.num_nodes, tri.num_nodes);
	s.setFromTriplets(tripletList.begin(), tripletList.end());

	return s;
}

SparseMatrix<double> P1FiniteElement::getFreeNodeStiffness() {
	return extractSparseMatrix(getStiffness(), fNodes, dof);
}

SparseMatrix<double> P1FiniteElement::getMass() {
	std::vector<T> tripletList;
	double value = 0;

	for (int t = 0; t < tri.num_edges; t++) {
		for (int i = 0; i < tri.edge_dim; i++) {
			for (int j = i; j < tri.edge_dim; j++) {
				if (i == j) {
					value = volumes(t) * 2 / (tri.edge_dim * (tri.edge_dim + 1));
				}
				else {
					value = volumes(t) / (tri.edge_dim * (tri.edge_dim + 1));
				}

				tripletList.push_back(T(tri.n4e(t, i), tri.n4e(t, j), value));
				if (i != j) tripletList.push_back(T(tri.n4e(t, j), tri.n4e(t, i), value));
			}
		}
	}
	SparseMatrix<double> m(tri.num_nodes, tri.num_nodes);
	m.setFromTriplets(tripletList.begin(), tripletList.end());

	return m;
}

SparseMatrix<double> P1FiniteElement::getFreeNodeMass() {
	return extractSparseMatrix(getMass(), fNodes, dof);
}

VectorXd P1FiniteElement::getRHSMassIntegrator(std::function<double(VectorXd&)> f) {
	VectorXd rhs = VectorXd::Zero(tri.num_nodes);
	for (int t = 0; t < tri.num_edges; t++) {
		for (int i = 0; i < tri.edge_dim; i++) {
			VectorXd v = bctr.col(t);
			rhs(tri.n4e(t, i)) += f(v) * volumes(t) / tri.edge_dim;
		}
	}
	return rhs;
}

VectorXd P1FiniteElement::getFreeNodeRHSMassIntegrator(std::function<double(VectorXd&)> f) {
	return extractDenseVector(getRHSMassIntegrator(f), fNodes, dof);
}

VectorXd P1FiniteElement::addDirichletNodes(VectorXd X) {
	return insertDenseVector(X, fNodes, tri.num_nodes);
}
