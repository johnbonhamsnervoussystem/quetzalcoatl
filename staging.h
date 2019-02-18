#include "common.h"
#include <Eigen/Dense>
#include "tei.h"
#include <vector>

#ifndef STAGING_H
#define STAGING_H

/*
  Routines for staging jobs are wrapped up here.  These should serve to 
  keep main clean and readable so that the parsing of jobs is not cluttered
  with generation of necessary elements.
*/
void molecular_hamiltonian( common& c) ;

void nnrep( common& com) ;

void nnrep( common& com, int& n, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> a ) ;

void pairing_hamiltonian( common& c) ;

#endif
