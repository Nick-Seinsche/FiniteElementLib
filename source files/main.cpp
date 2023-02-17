#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseLU>

#include "../header files/triangulation.h"
#include "../header files/P1FiniteElement.h"
#include "../header files/visual.h"
#include "../header files/gridfunction.h"

using Eigen::MatrixXd;
using namespace std;


// todo
// 3d mesh refine
// add neumann boundary support


void test() {
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    cout << m << endl;
}

int main(){ 
    cout << "calculating triangulation..." << endl;

    // create a triangulation object
    // by default this is the coarse unit square triangulation
    Triangulation tri;
    // do a couple of uniform refinements
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

    // calculate the stiffness matrix for the free nodes
    SparseMatrix<double> S = t.getFreeNodeStiffness();

    // calculate the right hand side stiffness matrix to incorperate non-zero
    // dirichlet boundary
    SparseMatrix<double> SR = extractVerticalSparseMatrix(t.getStiffness(), t.fNodes);

    // calculate the mass matrix on free nodes
    SparseMatrix<double> M = t.getFreeNodeMass();

    SparseMatrix<double> MR = extractVerticalSparseMatrix(t.getMass(), t.fNodes);

    // calculate the right hand side
    VectorXd B = t.getFreeNodeRHSMassIntegrator([](VectorXd x) { 
            return 10 * ((x[0] - 0.5) * (x[1] - 0.5));
        });

    // create dirichlet boundary vector
    VectorXd X_d = create([](VectorXd x) {
            return 0.01 * ((x[0] - 0.5) + (x[1] - 0.5));
        }, t.tri.c4n, t.tri.num_nodes);

    cout << "finite element ready!" << endl;
    cout << "solving linear system..." << endl;

    // Build -laplace(u) + 10 * u = f partial differential equation with 
    // non zero dirichlet boundary
    Eigen::SparseLU<SparseMatrix<double>> dec(S + 10 * M);
    //Eigen::PartialPivLU<MatrixXd> dec(S + 10 * M);
    VectorXd X_free = dec.solve(B - ((SR + 10 * MR) * X_d));

    cout << "solvied linear system!" << endl;

    // add the desired dirichlet boundary
    VectorXd X = t.addDirichletNodes(X_free) + X_d;

    cout << "visualizing..." << endl;

    tri.export_tri("mytriang");
    export_func(X, "mysol");

    trisurf();

    return 0;
}