#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#include "../../header files/P1FiniteElement.h"
#include "../../header files/visual.h"
#include "../../header files/io.h"

#include "../../header files/examples.h"

using namespace std;


void example_dist_heatflow() {
    double tau = 0.003;
    int T = 1;
    int N = int(T / tau);

    Matrix<double, 4, 2> c4n{ {0, 0}, {1, 0}, {1, 1}, {0, 1} };
    Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
    Matrix<int, 4, 2> Nb{ {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    Matrix<int, 0, 2> Db{  };

    Triangulation tri(c4n, n4e, Db, Nb);
    tri.uniform_refinement();
    //tri.uniform_refinement();
    //tri.uniform_refinement();
    //tri.uniform_refinement();

    P1FiniteElement t(tri);

    std::vector<int> ffNodes = t.fNodes;
    ffNodes.pop_back();
    SparseMatrix<double> S = extractSparseMatrix(t.getStiffness(), ffNodes, t.dof - 1);

    SparseMatrix<double> M = extractSparseMatrix(t.getMass(), ffNodes, t.dof - 1);

    Eigen::SparseLU<SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    
    SparseMatrix<double> I(t.dof - 1, t.dof - 1);
    I.setIdentity();

    //cout << I << endl;
    //cout << ((MatrixXd)S).determinant() << endl;
    solver.compute(S);
    SparseMatrix<double> SI = solver.solve(I);

    cout << S * SI << endl;

    VectorXd G = t.getNeumannBoundaryMassIntegrator([](VectorXd x) {
        return 0;
        });

    //VectorXd GR = G(t.fNodes);

    MatrixXd X(t.tri.num_nodes - 1, N);
    MatrixXd X_0(t.tri.num_nodes, N);

    VectorXd e = VectorXd::Zero(t.dof);
    e(t.dof - 1) = 1;

    X.col(0) = t.nodal_interpolant([](VectorXd x) {
        return (3 * (x[0] - 0.5) * (x[0] - 0.5) - 2 * (x[1] - 0.5)) / 2 - 0.3;
    }).block(0, 0, t.dof - 1, 1);

    X_0.block(0, 0, t.tri.num_nodes - 1, 1) = X.col(0);
    X_0(t.tri.num_nodes - 1, 0) = 0;

    X_0.col(0) = X_0.col(0) - (X_0.col(0).transpose() * t.m * VectorXd::Ones(t.dof)).value() * VectorXd::Ones(t.dof);
    X.col(0) = X_0.block(0, 0, t.tri.num_nodes - 1, 1);

    VectorXd RHS(t.tri.num_nodes);

    t.solver.compute(M * SI * M + tau * M);

    cout << VectorXd::Ones(t.dof) << endl;

    for (int k = 1; k < N; k++) {
        RHS = X.col(k - 1);
        X.col(k) = t.solver.solve(M * SI * M * RHS);
        X_0.block(0, k, t.tri.num_nodes - 1, 1) = X.col(k);
        X_0(t.tri.num_nodes - 1, k) = 0;
        //cout << (X_0.col(k).transpose() * t.m * VectorXd::Ones(t.dof)).value() << endl;
        //cout << (e.transpose() * t.m * VectorXd::Ones(t.dof)).value() << endl;
        //return;
        double qqq = (-(X_0.col(k).transpose() * t.m * VectorXd::Ones(t.dof)).value() / (e.transpose() * t.m * VectorXd::Ones(t.dof)).value());
        cout << "a " << -(X_0.col(k).transpose() * t.m * VectorXd::Ones(t.dof)).value() << endl;
        cout << "b " << (e.transpose() * t.m * VectorXd::Ones(t.dof)).value() << endl;
        //cout << (-(X_0.col(k).transpose() * t.m * VectorXd::Ones(t.dof)).value() / (e.transpose() * t.m * VectorXd::Ones(t.dof)).value()) << endl;
        X_0(t.tri.num_nodes - 1, k) = qqq;
        cout << "c " << X_0(t.tri.num_nodes - 1, k) << endl;
    }


    //X_0.block(0, 0, t.tri.num_nodes - 1, N) = X;
    //X_0.row(t.tri.num_nodes - 1) = MatrixXd::Zero(1, N);

    tri.export_tri("5_square_triang");
    IO::export_matrix_as_tsv(X_0, IO::SOLUTION, "5_dist_heat");

    //export_mat(X_0, "5_dist_heat");

    visual::trisurf("5_square_triang", "5_dist_heat");
    //visual::tricontour("mytriang", "mysol");
}