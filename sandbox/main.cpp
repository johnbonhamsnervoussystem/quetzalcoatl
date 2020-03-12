#include <Eigen/Dense>
#include "hamiltonian.h"

int main(void){
  hamiltonian::molecular t;
  Eigen::VectorXd a;
  Eigen::MatrixXd c;
  a.resize(2);
  c.resize(2, 3);
  a(0) = 1.0;
  a(1) = 1.0;
//
  c(0, 0) = 0.0;
  c(0, 1) = 0.0;
  c(0, 2) = 0.0;
 //
  c(1, 0) = 0.0;
  c(1, 1) = 0.0;
  c(1, 2) = 0.74;

//  t = hamiltonian::molecular(a, c);
  return 0;
  }
