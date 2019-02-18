#include "constants.h"
#include <Eigen/Core>
#include "evalm.h"
#include "nbodyint.h"
#include <iostream>
#include "qtzio.h"
#include "pairing.h"
#include "util.h"

template <class matrix>
void pairing<matrix>::contract( matrix& m){

  if ( itype == 1){
    G = -U*m ;
    }

  return ;

} ;

template void pairing<Eigen::MatrixXd>::contract( Eigen::MatrixXd&) ;

template void pairing<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&) ;

template <class matrix>
void pairing<matrix>::contract( matrix& m, matrix& n){

  return ;

} ;

template void pairing<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;

template void pairing<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;

