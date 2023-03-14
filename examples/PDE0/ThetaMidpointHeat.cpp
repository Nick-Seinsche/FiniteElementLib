#include <iostream>
#include <Eigen/Dense>

#include "../../header files/P1FiniteElement.h"
#include "../../header files/visual.h"
#include "../../header files/io.h"

#include "../../header files/examples.h"

using namespace std;


void exampleThetaMidpointHeat() {
    double tau = 0.004;
    int T = 1;
    int N = int(T / tau);
    double kappa = 0.3;
    float theta = 0.5;

    Matrix<double, 4, 2> c4n{ {0, 0}, {1, 0}, {1, 1}, {0, 1} };
    Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
    Matrix<int, 3, 2> Db{ {1, 2}, {2, 3}, {3, 0} };
    Matrix<int, 1, 2> Nb{ {0, 1} };

    Triangulation tri(c4n, n4e, Db, Nb);
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();

    P1FiniteElement t(tri);

    SparseMatrix<double> S = t.getStiffness();
    SparseMatrix<double> M = t.getMass();

    VectorXd G = t.getNeumannBoundaryMassIntegrator([](VectorXd x) {
        return 2;
        });

    VectorXd GR = G(t.fNodes);

    t.solver.compute(M + theta * kappa * tau * S);

    MatrixXd X(t.tri.num_nodes, N);

    X.col(0) = t.nodal_interpolant([](VectorXd x) {
        return 0;
        });

    VectorXd RHS(t.tri.num_nodes);

    for (int k = 1; k < N; k++) {
        RHS = X.col(k - 1);
        X.col(k) = insertDenseVector(t.solver.solve(M * RHS(t.fNodes) + kappa * tau * GR - kappa * tau * (1 - theta) * S * RHS(t.fNodes)),
            t.fNodes, t.tri.num_nodes, 0);
    }

    tri.export_tri("1_square_triang");
    //export_mat(X, "1_theta_heat");
    IO::export_matrix_as_tsv(X, IO::SOLUTION, "1_theta_heat");

    visual::trisurf("1_square_triang", "1_theta_heat");
    //visual::tricontour("mytriang", "mysol");
}