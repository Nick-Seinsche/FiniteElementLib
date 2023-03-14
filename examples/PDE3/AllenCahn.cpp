#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>

#include "../../header files/P1FiniteElement.h"
#include "../../header files/visual.h"
#include "../../header files/io.h"

#include "../../header files/examples.h"


VectorXd calculateWp(P1FiniteElement& p, VectorXd X) {
    auto Wp = [](double x) {return 4 * x * (x * x - 1); };
    VectorXd rhs = VectorXd::Zero(p.tri.num_nodes);
    VectorXd v;
    for (int t = 0; t < p.tri.num_edges; t++) {
        for (int i = 0; i < p.tri.edge_dim; i++) {
            v = p.bctr.col(t);
            rhs(p.tri.n4e(t, i)) += Wp((X(p.tri.n4e(t,0)) + X(p.tri.n4e(t, 1)) + X(p.tri.n4e(t, 2)))/3) * p.volumes(t) / p.tri.edge_dim;
        }
    }
    return rhs;
}

SparseMatrix<double> calculateWpp(P1FiniteElement& p, VectorXd X) {
    auto Wpp = [](double x) {return 12 * x * x - 4; };
    std::vector<T> tripletList;
    double value = 0;

    for (int t = 0; t < p.tri.num_edges; t++) {
        for (int i = 0; i < p.tri.edge_dim; i++) {
            for (int j = i; j < p.tri.edge_dim; j++) {
                if (i == j) {
                    value = Wpp((X(p.tri.n4e(t, 0)) + X(p.tri.n4e(t, 1)) + X(p.tri.n4e(t, 2))) / 3) * p.volumes(t) * 2 / (p.tri.edge_dim * (p.tri.edge_dim + 1));
                }
                else {
                    value = Wpp((X(p.tri.n4e(t, 0)) + X(p.tri.n4e(t, 1)) + X(p.tri.n4e(t, 2))) / 3) * p.volumes(t) / (p.tri.edge_dim * (p.tri.edge_dim + 1));
                }

                tripletList.push_back(T(p.tri.n4e(t, i), p.tri.n4e(t, j), value));
                if (i != j) tripletList.push_back(T(p.tri.n4e(t, j), p.tri.n4e(t, i), value));
            }
        }
    }
    SparseMatrix<double> w(p.tri.num_nodes, p.tri.num_nodes);
    w.setFromTriplets(tripletList.begin(), tripletList.end());

    return w;
}


void example_allen_cahn() {
    double tau = 0.001;
    double Time = 0.12;
    int N = int(Time / tau);
    double eps = 3;

    Matrix<double, 4, 2> c4n{ {0, 0}, {1, 0}, {1, 1}, {0, 1} };
    Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
    Matrix<int, 4, 2> Db{ {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    Matrix<int, 0, 2> Nb{  };

    Triangulation tri(c4n, n4e, Db, Nb);
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();

    P1FiniteElement t(tri);

    SparseMatrix<double> S = t.getStiffness();
    SparseMatrix<double> M = t.getMass();

    MatrixXd X(t.tri.num_nodes, N);
    VectorXd X_old(t.tri.num_nodes);

    VectorXd U_old = VectorXd::Zero(t.dof);
    VectorXd U(t.dof);
    VectorXd rhs(t.dof);

    SparseMatrix<double> WppMat(t.tri.num_nodes, t.tri.num_nodes);
    SparseMatrix<double> WppMatf(t.dof, t.dof);
    VectorXd WpVec(t.tri.num_nodes);

    X.col(0) = t.nodal_interpolant([](VectorXd x) {
        if (abs(x[0] - 0.5) + abs(x[1] - 0.5) <= 0.27) {
            return 1;
        }
        else {
            return -1;
        }
    });

    for (int k = 1; k < N; k++) {
        std::cout << (float) k / N << std::endl;
        X_old = X.col(k - 1);
        U = X_old(t.fNodes);

        do {
            U_old = U;

            WppMat = calculateWpp(t, insertDenseVector(U_old, t.fNodes, t.tri.num_nodes, 0));

            WppMatf = extractSparseMatrix(WppMat, t.fNodes, t.dof);
            WpVec = calculateWp(t, insertDenseVector(U_old, t.fNodes, t.tri.num_nodes, 0));

            rhs = S * U_old + WpVec(t.fNodes) / (eps * eps) + M * (U_old - X_old(t.fNodes)) / tau;

            Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
            cg.compute(S + M / tau + WppMatf / (eps * eps));
            U = cg.solve((S + M / tau + WppMatf / (eps * eps)) * U_old - rhs);

            //std::cout << (U - U_old).norm() << std::endl;
        } while ((U - U_old).norm() > 0.001);

        //std::cout << "newton" << std::endl;

        X.col(k) = insertDenseVector(U, t.fNodes, t.tri.num_nodes, 0);
    }

    tri.export_tri("3_square_triang");
    //export_mat(X, "3_allen_cahn");
    IO::export_matrix_as_tsv(X, IO::SOLUTION, "3_allen_cahn");
    
    //visual::tricontour("3_square_triang", "3_allen_cahn");
    visual::trisurf("3_square_triang", "3_allen_cahn");
}