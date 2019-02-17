#include "constants.h"
#include <Eigen/Core>
#include "evalm.h"
#include "nbodyint.h"
#include <iostream>
#include "r12.h"
#include "util.h"

template <class matrix>
void r12<matrix>::contract( matrix& m){

  if ( itype == 1){

/*
  rhf
*/
    transform( 2, x, m) ;
    ctr2er( *ints, m, G, dim) ;
    transform( 2, xi, m) ;
  } else if ( itype == 3) {

/*
  ghf
*/
    transform( 2, x, m) ;
    ctr2eg( *ints, m, G, dim/2) ;
    transform( 2, xi, m) ;
    }

  return ;

} ;

template void r12<Eigen::MatrixXd>::contract( Eigen::MatrixXd&) ;

template void r12<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&) ;

template <class matrix>
void r12<matrix>::contract( matrix& m, matrix& n){

  if ( itype == 2){

/*
  uhf
*/
    transform( 2, x, m) ;
    transform( 2, x, n) ;
    ctr2eu( *ints, m, n, G, dim) ;
    transform( 2, xi, m) ;
    transform( 2, xi, n) ;

  } else if ( itype == 4){

/*
  rhfb
*/
    transform( 2, x, m) ;
    transform( 2, x, n) ;
    ctr2er( *ints, m, G.block( 0, 0, dim/2, dim/2), dim/2) ;
    ctrPairr( *ints, n, G.block( 0, dim/2, dim/2, dim/2), dim/2) ;
    transform( 2, xi, m) ;
    transform( 2, xi, n) ;
    transform( 2, x, G.block( 0, 0, dim/2, dim/2)) ;
    transform( 2, x, G.block( 0, dim/2, dim/2, dim/2)) ;
    G.block( 0, dim/2, dim/2, dim/2) /= d2 ;
    G.block( dim/2, dim/2, dim/2, dim/2) = -G.block( 0, 0, dim/2, dim/2) ;
    G.block( dim/2, 0, dim/2, dim/2) = -G.block( 0, dim/2, dim/2, dim/2) ;
 
  } else if ( itype == 6){

/*
  ghfb
*/
    transform( 2, x, m) ;
    transform( 2, x, n) ;
    ctr2eg( *ints, m, G.block( 0, 0, dim/2, dim/2), dim/2) ;
    ctrPairg( *ints, n, G.block( 0, dim/2, dim/2, dim/2), dim/2) ;
    transform( 2, xi, m) ;
    transform( 2, xi, n) ;
    transform( 2, x, G.block( 0, 0, dim/2, dim/2)) ;
    transform( 2, x, G.block( 0, dim/2, dim/2, dim/2)) ;
    G.block( 0, dim/2, dim/2, dim/2) /= d2 ;
    G.block( dim/2, dim/2, dim/2, dim/2) = -G.block( 0, 0, dim/2, dim/2) ;
    G.block( dim/2, 0, dim/2, dim/2) = -G.block( 0, dim/2, dim/2, dim/2) ;
 
    }

  return ;

} ;

template void r12<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;

template void r12<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;

