#include "common.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "nbodyint.h"
#include "tei.h"

#ifndef HFROUT_H
#define HFROUT_H

/*
class SlaterDetDriver {
  * Driver class for solving for SlaterDeterminants.
  * If no special case is set choose the most general version.
  * We need the intergral engines/atoms or matrices with the values
  * We need algorithm options
  private:
    nalpha
    nbeta
    nbasis
    convergence

    
  public:
    
  
}


class BogolyubovDriver {
  * Driver class for solving for Hartree-Fock Bogolyuboc wavefunctions.
  * If no special case is set choose the most general version.
  private:
  public:
  
}


class RealResHfDia {
  private:
    nbasis
    density
    fock
    nbody
    hcore
    convergence
    energy
    eigensolver
    iteration_limit
  public:
    void execute(void) {
       * Perform the minimization
    }

  }
*/

/*
double rhfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nele, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, bool ls, double lvs, const int& maxit, const double& thresh) ;
#endif
