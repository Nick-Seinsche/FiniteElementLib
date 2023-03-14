#include <iostream>
#include <cmath>
#include <Eigen/Dense>

#include "../../header files/P1FiniteElement.h"
#include "../../header files/io.h"
#include "../../header files/visual.h"

#include "../../header files/examples.h"

using namespace std;

void examplePoisson() {
    cout << "calculating triangulation..." << endl;

    Matrix<double, 8, 2> c4n{ {0.0, 0.0}, {0.5, 0}, {1, 0}, {0.8, 0.5}, {1, 1}, {0.5, 1}, {0, 1}, {0.2, 0.5} };
    Matrix<int, 6, 3> n4e{ { 0, 1, 7 }, { 1, 2, 3 }, {1, 3, 7}, {3, 4, 5}, {7, 3, 5}, {7, 5, 6} };
    Matrix<int, 4, 2> Nb{ {0, 7}, {6, 7}, {2, 3}, {3, 4} };
    Matrix<int, 4, 2> Db{ {0, 1}, {1, 2}, {6, 5}, {5, 4} };

    //Matrix<double, 4, 2> c4n{ {0.4, 0.4}, {1, 0}, {1, 1}, {0, 1} };
    //Matrix<int, 2, 3> n4e{ { 0, 1, 2 }, { 0, 2, 3 } };
    //Matrix<int, 2, 2> Db{ {0, 1}, {2, 3} };
    //Matrix<int, 2, 2> Nb{ {1, 2}, {3, 0} };

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
    //tri.uniform_refinement();
    //tri.uniform_refinement();

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

    // calculate the right hand side
    VectorXd B = t.getRHSMassIntegrator([](VectorXd x) {
        return 0;
        });
    VectorXd BR = B(t.fNodes);

    VectorXd G = t.getNeumannBoundaryMassIntegrator([](VectorXd x) {
        return 0 ;
        });

    VectorXd GR = G(t.fNodes);

    // create dirichlet boundary vector
    VectorXd X_d = t.nodal_interpolant([](VectorXd x) {
        return sin(3.14159 * x[0]);
    });

    cout << "finite element ready!" << endl;
    cout << "solving linear system..." << endl;

    // Build -laplace(u) + u = f partial differential equation with 
    // non zero dirichlet boundary

    t.solver.compute(S);

    VectorXd X_free = t.solver.solve(BR + GR - SR * X_d);

    cout << "solvied linear system!" << endl;

    // add the desired dirichlet boundary
    VectorXd X = insertDenseVector(X_free, t.fNodes, tri.num_nodes, 0) + X_d;

    cout << "visualizing..." << endl;

    tri.export_tri("00_hex_triang");
    IO::export_vector_as_tsv(X, IO::SOLUTION, "00_poisson");

    visual::trisurf("00_hex_triang", "00_poisson");
    //visual::tricontour("mytriang", "mysol");

}