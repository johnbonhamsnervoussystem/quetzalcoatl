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

  if ( itype == 4){
    /*
      rhfb 
    */
    G.setZero() ;
    G.block( 0, 0, dim, dim) = -U*m ;
    G.block( dim, dim, dim, dim) = -G.block( 0, 0, dim, dim).conjugate() ;
    G.topRightCorner( dim, dim).setIdentity() ;
    G.topRightCorner( dim, dim) *= -U*n.trace() ;
    G.block( dim, 0, dim, dim) = -G.block( 0, dim, dim, dim).conjugate() ;
    }

  return ;

} ;

template void pairing<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;

template void pairing<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;

