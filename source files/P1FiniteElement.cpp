#include <iostream>
#include <thread>
#include <mutex>

#include "../header files/P1FiniteElement.h"

//std::mutex mu;

int __sfactorial(int k) {
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

SparseMatrix<double> extractVerticalSparseMatrix(SparseMatrix<double> M, std::vector<int> nodes, int dof) {
	SparseMatrix<double> M_red(dof, M.cols());
	for (int k = 0; k < M.outerSize(); k++)
		for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
			auto rw = std::find(nodes.begin(), nodes.end(), it.row());
			if (rw != nodes.end()) {
				M_red.insert(rw - nodes.begin(), it.col()) = it.value();
			}
		}
	return M_red;
}

VectorXd insertDenseVector(VectorXd V, std::vector<int> fNodes, int full, double value) {
	VectorXd V_full(full);
	int j = 0;
	for (int i = 0; i < full; i++) {
		if (std::find(fNodes.begin(), fNodes.end(), i) != fNodes.end()) {
			V_full(i) = V(j);
			j++;
		} else {
			V_full(i) = value;
		}
	}
	return V_full;
}

void P1FiniteElement::InitP1FiniteElementThread(Triangulation& t, int index, int num_threads) {
	MatrixXd X_T(tri.node_dim + 1, tri.node_dim + 1);
	MatrixXd I_d(tri.node_dim + 1, tri.node_dim);
	I_d.row(0) = MatrixXd::Zero(1, tri.node_dim);
	I_d.block(1, 0, tri.node_dim, tri.node_dim)
		= MatrixXd::Identity(tri.node_dim, tri.node_dim);

	VectorXd bctri;

	for (int i = 0; i < tri.num_edges; i++) {
		if (i % num_threads != index) continue;
		bctri = VectorXd::Zero(tri.node_dim);
		for (int j = 0; j < tri.edge_dim; j++) {
			X_T(0, j) = 1;
			X_T.block(1, j, tri.node_dim, 1) = tri.c4n.row(tri.n4e(i, j)).transpose();

			bctri += tri.c4n.row(tri.n4e(i, j)).transpose();
		}

		volumes(i) = X_T.determinant() / __sfactorial(tri.node_dim);
		grads[i] = (X_T.inverse() * I_d).transpose();

		bctr.col(i) = bctri / tri.edge_dim;
	}
}

P1FiniteElement::P1FiniteElement(Triangulation tri) {
	this->tri = tri;
	this->grads = new MatrixXd[tri.num_edges];
	this->volumes = VectorXd(tri.num_edges);
	this->bctr = MatrixXd(tri.node_dim, tri.num_edges);

	std::set<int> sdb;
	for (auto k : tri.Db.reshaped()) {
		sdb.insert(k);
	}
	std::vector<int> fNodes;
	int j = 0;
	for (int i = 0; i < tri.num_nodes; i++) {
		if (sdb.find(i) == sdb.end()) {
			this->ifNodes.push_back(j);
			j++;
			fNodes.push_back(i);
			this->dof++;
		}
		else {
			this->ifNodes.push_back(-1);
		}
	}
	this->fNodes = fNodes;

	VectorXd bctrnbi;
	VectorXd vol_Nb(tri.Nb.rows());
	MatrixXd bctr_Nb(tri.node_dim, tri.Nb.rows());

	for (int i = 0; i < tri.Nb.rows(); i++) {
		bctrnbi = VectorXd::Zero(tri.node_dim);
		for (int j = 0; j < tri.edge_dim - 1; j++) {
			bctrnbi += tri.c4n.row(tri.Nb(i, j)).transpose();
		}
		bctr_Nb.col(i) = bctrnbi / (tri.edge_dim - 1);

		if (tri.node_dim == 1) vol_Nb(i) = 1;
		else if (tri.node_dim == 2) vol_Nb(i) = (tri.c4n.row(tri.Nb(i, 1)) - tri.c4n.row(tri.Nb(i, 0))).stableNorm();
		//else if (tri.node_dim == 3) {
		//	temp = tri.c4n.row(tri.Nb(i, 2)) - tri.c4n.row(tri.Nb(i, 0));
		//	vol_Nb(i) = temp.cross(tri.c4n.row(tri.Nb(i, 1)) - tri.c4n.row(tri.Nb(i, 0))).stableNorm() / 2;
		//}
		else std::cerr << "volumes not implemented for d = 3 yet." << std::endl;
	}

	this->volumes_Nb = vol_Nb;
	this->bctr_Nb = bctr_Nb;

	int num_threads = 16;
	std::vector<std::thread> threads;

	for (int i = 0; i < num_threads; i++) {
		threads.push_back(std::thread(&P1FiniteElement::InitP1FiniteElementThread, this, std::ref(tri), i, num_threads));
	}
	for (std::thread& v : threads) v.join();
}

void P1FiniteElement::getStiffnessThread(SparseMatrix<double>& s, SparseMatrix<double>& sff,
	SparseMatrix<double>& sr, int index, int num_threads) {
	
	std::vector<T> tripletList_s;
	std::vector<T> tripletList_sff;
	std::vector<T> tripletList_sr;
	double value = 0;
	for (int t = 0; t < tri.num_edges; t++) {
		if (t % num_threads != index) continue;
		for (int i = 0; i < tri.edge_dim; i++) {
			for (int j = i; j < tri.edge_dim; j++) {
				value = volumes(t) * (grads[t].col(i).transpose() * grads[t].col(j)).value();
				
				tripletList_s.push_back(T(tri.n4e(t, i), tri.n4e(t, j), value));
				if (i != j) tripletList_s.push_back(T(tri.n4e(t, j), tri.n4e(t, i), value));

				if (ifNodes[tri.n4e(t, i)] != -1 && ifNodes[tri.n4e(t, j)] != -1) {
					tripletList_sff.push_back(T(ifNodes[tri.n4e(t, i)], ifNodes[tri.n4e(t, j)], value));
					if (i != j) tripletList_sff.push_back(T(ifNodes[tri.n4e(t, j)], ifNodes[tri.n4e(t, i)], value));
				}

				if (ifNodes[tri.n4e(t, i)] != -1) {
					tripletList_sr.push_back(T(ifNodes[tri.n4e(t, i)], tri.n4e(t, j), value));
				}

				if (ifNodes[tri.n4e(t, j)] != -1 && i != j) {
					tripletList_sr.push_back(T(ifNodes[tri.n4e(t, j)], tri.n4e(t, i), value));
				}
			}
		}
	}
	s = SparseMatrix<double>(tri.num_nodes, tri.num_nodes);
	s.setFromTriplets(tripletList_s.begin(), tripletList_s.end());

	sff = SparseMatrix<double>(dof, dof);
	sff.setFromTriplets(tripletList_sff.begin(), tripletList_sff.end());

	sr = SparseMatrix<double>(dof, tri.num_nodes);
	sr.setFromTriplets(tripletList_sr.begin(), tripletList_sr.end());
}


SparseMatrix<double> P1FiniteElement::getStiffness() {
	int num_threads = 16;
	
	this->s = SparseMatrix<double>(tri.num_nodes, tri.num_nodes);
	this->sff = SparseMatrix<double>(dof, dof);
	this->sr = SparseMatrix<double>(dof, tri.num_nodes);
	
	std::vector<SparseMatrix<double>> s_th(num_threads);
	std::vector<SparseMatrix<double>> sff_th(num_threads);
	std::vector<SparseMatrix<double>> sr_th(num_threads);
	
	std::vector<std::thread> threads;
	
	for (int i = 0; i < num_threads; i++) {
		threads.push_back(std::thread(&P1FiniteElement::getStiffnessThread, this, std::ref(s_th[i]), std::ref(sff_th[i]), 
			std::ref(sr_th[i]), i, num_threads));
	}

	for (std::thread& v : threads) v.join();
	
	for (int i = 0; i < num_threads; i++) {
		this->s += s_th[i];
		this->sff += sff_th[i];
		this->sr += sr_th[i];
	}

	return sff;
}

void P1FiniteElement::getMassThread(SparseMatrix<double>& m, SparseMatrix<double>& mff,
	SparseMatrix<double>& mr, int index, int num_threads) {
	
	std::vector<T> tripletList_m;
	std::vector<T> tripletList_mff;
	std::vector<T> tripletList_mr;
	double value = 0;
	for (int t = 0; t < tri.num_edges; t++) {
		if (t % num_threads != index) continue;
		for (int i = 0; i < tri.edge_dim; i++) {
			for (int j = i; j < tri.edge_dim; j++) {
				if (i == j) {
					value = volumes(t) * 2 / (tri.edge_dim * (tri.edge_dim + 1));
				}
				else {
					value = volumes(t) / (tri.edge_dim * (tri.edge_dim + 1));
				}

				tripletList_m.push_back(T(tri.n4e(t, i), tri.n4e(t, j), value));
				if (i != j) tripletList_m.push_back(T(tri.n4e(t, j), tri.n4e(t, i), value));

				if (ifNodes[tri.n4e(t, i)] != -1 && ifNodes[tri.n4e(t, j)] != -1) {
					tripletList_mff.push_back(T(ifNodes[tri.n4e(t, i)], ifNodes[tri.n4e(t, j)], value));
					if (i != j) tripletList_mff.push_back(T(ifNodes[tri.n4e(t, j)], ifNodes[tri.n4e(t, i)], value));
				}

				if (ifNodes[tri.n4e(t, i)] != -1) {
					tripletList_mr.push_back(T(ifNodes[tri.n4e(t, i)], tri.n4e(t, j), value));
				}

				if (ifNodes[tri.n4e(t, j)] != -1 && i != j) {
					tripletList_mr.push_back(T(ifNodes[tri.n4e(t, j)], tri.n4e(t, i), value));
				}
			}
		}
	}
	m = SparseMatrix<double>(tri.num_nodes, tri.num_nodes);
	m.setFromTriplets(tripletList_m.begin(), tripletList_m.end());

	mff = SparseMatrix<double>(dof, dof);
	mff.setFromTriplets(tripletList_mff.begin(), tripletList_mff.end());

	mr = SparseMatrix<double>(dof, tri.num_nodes);
	mr.setFromTriplets(tripletList_mr.begin(), tripletList_mr.end());
}

SparseMatrix<double> P1FiniteElement::getMass() {
	int num_threads = 16;

	this->m = SparseMatrix<double>(tri.num_nodes, tri.num_nodes);
	this->mff = SparseMatrix<double>(dof, dof);
	this->mr = SparseMatrix<double>(dof, tri.num_nodes);

	std::vector<SparseMatrix<double>> m_th(num_threads);
	std::vector<SparseMatrix<double>> mff_th(num_threads);
	std::vector<SparseMatrix<double>> mr_th(num_threads);

	std::vector<std::thread> threads;

	for (int i = 0; i < num_threads; i++) {
		threads.push_back(std::thread(&P1FiniteElement::getMassThread, this, std::ref(m_th[i]), std::ref(mff_th[i]),
			std::ref(mr_th[i]), i, num_threads));
	}
	for (std::thread& v : threads) v.join();

	for (int i = 0; i < num_threads; i++) {
		this->m += m_th[i];
		this->mff += mff_th[i];
		this->mr += mr_th[i];
	}

	return mff;
}

VectorXd P1FiniteElement::getRHSMassIntegrator(std::function<double(VectorXd&)> f) {
	VectorXd rhs = VectorXd::Zero(tri.num_nodes);
	VectorXd v;
	for (int t = 0; t < tri.num_edges; t++) {
		for (int i = 0; i < tri.edge_dim; i++) {
			v = bctr.col(t);
			rhs(tri.n4e(t, i)) += f(v) * volumes(t) / tri.edge_dim;
		}
	}
	return rhs;
}

VectorXd P1FiniteElement::getNeumannBoundaryMassIntegrator(std::function<double(VectorXd&)> g) {
	VectorXd rhs = VectorXd::Zero(tri.num_nodes);
	for (int t = 0; t < tri.Nb.rows(); t++) {
		for (int i = 0; i < tri.edge_dim - 1; i++) {
			VectorXd v = bctr_Nb.col(t);
			rhs(tri.Nb(t, i)) += g(v) * volumes_Nb(t) / (tri.edge_dim - 1);
		}
	}
	return rhs;
}

VectorXd P1FiniteElement::nodal_interpolant(std::function<double(VectorXd)> f) {
	VectorXd gridfunction(tri.num_nodes);
	for (int i = 0; i < tri.num_nodes; i++) {
		gridfunction(i) = f(tri.c4n.row(i).transpose());
	}
	return gridfunction;
}