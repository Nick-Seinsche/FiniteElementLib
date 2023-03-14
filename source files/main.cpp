#include <iostream>
#include <Eigen/Dense>

#include "../header files/examples.h"

#include <thread>

using Eigen::MatrixXd;
using namespace std;

// todo
// 3d mesh refine / 1d mesh refine
// gradient flow
// restructuring of lib
// faster extract insert !!


void test() {
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    cout << m << endl;
}

int main(){ 
    //examplePoisson();
    //exampleEigenvalue();
   //exampleThetaMidpointHeat();
    //example_dist_heatflow();
    //visual::tricontour("3_square_triang", "3_allen_cahn");
    example_allen_cahn();
    //example_cahn_hilliard();
    return 0;
}

