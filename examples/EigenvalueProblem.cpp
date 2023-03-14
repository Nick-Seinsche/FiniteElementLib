#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseLU>

#include "../header files/P1FiniteElement.h"
#include "../header files/visual.h"
#include "../header files/io.h"

#include "../header files/examples.h"

using namespace std;

void exampleEigenvalue() {
    cout << "calculating triangulation..." << endl;

    Matrix<double, 4, 2> c4n{ {0, 0}, {1, 0}, {1, 1}, {0, 1} };
    Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
    Matrix<int, 2, 2> Db{ {0, 1}, {2, 3} };
    Matrix<int, 2, 2> Nb{ {1, 2}, {3, 0} };

    // create a triangulation object
    // by default this is the coarse unit square triangulation
    Triangulation tri(c4n, n4e, Db, Nb);

    // do a couple of uniform refinements
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();
    tri.uniform_refinement();

    cout << "Triangulation ready!" << endl;
    cout << "Setting up finite element..." << endl;

    // create the P1FiniteElement Class by passing the triangulation
    P1FiniteElement t(tri);

    cout << "degress of freedom: " << t.dof << endl;

    // todo: parallelization
    // calculate the stiffness matrix for the free nodes
    SparseMatrix<double> S = t.getStiffness();

    cout << "calculated stiffness" << endl;

    // calculate the right hand side stiffness matrix to incorperate non-zero
    // dirichlet boundary
    SparseMatrix<double> SR = t.sr;

    // calculate the mass matrix on free nodes
    SparseMatrix<double> M = t.getMass();

    cout << "calculated mass" << endl;

    SparseMatrix<double> MR = t.mr;

    // calculate the right hand side
    VectorXd B = t.getRHSMassIntegrator([](VectorXd x) {
        return 0;
        });
    VectorXd BR = B(t.fNodes);

    VectorXd G = t.getNeumannBoundaryMassIntegrator([](VectorXd x) {
        return 0.5 * (x[0] - 0.5);
        });

    VectorXd GR = G(t.fNodes);

    // create dirichlet boundary vector
    VectorXd X_d = t.nodal_interpolant([](VectorXd x) {
        return 0.5 * (1 - x[0]);
    });

    cout << "finite element ready!" << endl;
    cout << "solving linear system..." << endl;

    // Build -laplace(u) + u = f partial differential equation with 
    // non zero dirichlet boundary

    t.solver.compute(S + 10 * M);

    VectorXd X_free = t.solver.solve(BR + GR - ((SR + 10 * MR) * X_d));

    cout << "solvied linear system!" << endl;

    // add the desired dirichlet boundary
    VectorXd X = insertDenseVector(X_free, t.fNodes, tri.num_nodes, 0) + X_d;

    cout << "visualizing..." << endl;

    tri.export_tri("0_square_triang");
    //export_func(X, "0_manta_ray");
    IO::export_vector_as_tsv(X, IO::SOLUTION, "0_manta_ray");

    visual::trisurf("0_square_triang", "0_manta_ray");
    //visual::tricontour("mytriang", "mysol");
    
}