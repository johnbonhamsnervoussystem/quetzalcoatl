#include "constants.h"
#include <Eigen/Core>
#include "evalm.h"
#include "nbodyint.h"
#include <iostream>
#include "qtzio.h"
#include "r12.h"
#include "util.h"

template <class matrix>
void r12<matrix>::contract( matrix& m){

  if ( itype == 1){

/*
  rhf
*/
    transform( 1, x, m) ;
    ctr2er( *ints, m, G, dim) ;
    transform( 1, xi, m) ;
    transform( 0, x, G) ;
  } else if ( itype == 3) {

/*
  ghf
*/
    transform( 1, x, m) ;
    ctr2eg( *ints, m, G, dim) ;
    transform( 1, xi, m) ;
    transform( 0, x, G) ;
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
    transform( 1, x, m) ;
    transform( 1, x, n) ;
    ctr2eu( *ints, m, n, G, dim) ;
    transform( 1, xi, m) ;
    transform( 1, xi, n) ;
    transform( 0, x, G) ;

  } else if ( itype == 4){

/*
  rhfb
*/
    transform( 1, x, m) ;
    transform( 1, x, n) ;
    ctr2er( *ints, m, G.block( 0, 0, dim, dim), dim) ;
    ctrPairr( *ints, n, G.block( 0, dim, dim, dim), dim) ;
    transform( 1, xi, m) ;
    transform( 1, xi, n) ;
    transform( 0, x, G.block( 0, 0, dim, dim)) ;
    transform( 0, x, G.block( 0, dim, dim, dim)) ;
    G.block( 0, dim, dim, dim) /= two ;
    G.block( dim, dim, dim, dim) = -G.block( 0, 0, dim, dim).conjugate() ;
    G.block( dim, 0, dim, dim) = -G.block( 0, dim, dim, dim).conjugate() ;
 
  } else if ( itype == 6){

/*
  ghfb
*/
    transform( 1, x, m) ;
    transform( 1, x, n) ;
    ctr2eg( *ints, m, G.block( 0, 0, 2*dim, 2*dim), dim) ;
    ctrPairg( *ints, n, G.block( 0, 2*dim, 2*dim, 2*dim), dim) ;
    transform( 1, xi, m) ;
    transform( 1, xi, n) ;
    transform( 0, x, G.block( 0, 0, 2*dim, 2*dim)) ;
    transform( 0, x, G.block( 0, 2*dim, 2*dim, 2*dim)) ;
    G.block( 0, 2*dim, 2*dim, 2*dim) /= two ;
    G.block( 2*dim, 2*dim, 2*dim, 2*dim) = -G.block( 0, 0, 2*dim, 2*dim).conjugate() ;
    G.block( 2*dim, 0, 2*dim, 2*dim) = -G.block( 0, 2*dim, 2*dim, 2*dim).conjugate() ;
 
    }

  return ;

} ;

template void r12<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;

template void r12<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;

