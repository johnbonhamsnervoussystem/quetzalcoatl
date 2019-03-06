#include "constants.h"
#include <Eigen/Core>
#include "evalm.h"
#include "nbodyint.h"
#include <iostream>
#include "qtzio.h"
#include "qtzcntrl.h"
#include "hubbard.h"
#include "util.h"

template <class matrix>
void hubbard<matrix>::contract( matrix& m) {
  int i ;

  G.setZero() ;
  if ( itype == 1) {
    for( i=0; i < dim; i++) {
      G( i, i) = U*m( i, i) ;
      }
  } else if ( itype == 3){
    /*
      ghfb
      Gaa | Gab      pbb | -pab
      --------- = U* ----------
      Gba | Gbb      -pba| paa
    */
    G.setZero() ;
    for( i=0; i < dim; i++) {
      G( i, i) = U*m( i, i) ;
      G( dim + i, dim + i) = U*m( dim + i, dim + i) ;
      G( i + dim, i) = -U*m( dim + i, i) ;
      G( i, dim + i) = -U*m( i, dim + i) ;
      }
  } else {
    qtzcntrl::shutdown( " Unrecognized type in Hubbard::contract") ;
    }

  return ;

} ;

template void hubbard<Eigen::MatrixXd>::contract( Eigen::MatrixXd&) ;

template void hubbard<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&) ;

template <class matrix>
void hubbard<matrix>::contract( matrix& m, matrix& n){
  int i ;

  if ( itype == 4){
    /*
      rhfb
      K^{ab} = -K^{ba}
    */
    G.setZero() ;
    for( i=0; i < dim; i++) {
      G( i, i) = U*m( i, i) ;
      G( i, dim+ i) = -U*n( i, i) ;
      }
    G.block( dim, dim, dim, dim) = -G.block( 0, 0, dim, dim).conjugate() ;
    G.block( dim, 0, dim, dim) = -G.block( 0, dim, dim, dim).conjugate() ;
  } else if ( itype == 6){
    /*
      ghfb
      Gaa | Gab      pbb | -pba
      --------- = U* ----------
      Gba | Gbb      -pab| paa

      Daa | Dab       0  | kba
      --------- = U* ----------
      Dba | Dbb       kab|  0 
    */
    G.setZero() ;
    for( i = 0; i < dim; i++) {
      G( i, i) = U*m( dim + i, dim + i) ;
      G( dim + i, dim + i) = U*m( i, i) ;
      G( i + dim, i) = -U*m( i, dim + i) ;
      G( i, dim + i) = -U*m( dim + i, i) ;
      G( i, 3*dim + i) = U*n( dim + i, i) ;
      G( dim + i, 2*dim + i) = U*n( i, dim + i) ;
      }
    G.block( 2*dim, 2*dim, 2*dim, 2*dim) = -G.block( 0, 0, 2*dim, 2*dim).conjugate() ;
    G.block( 2*dim, 0, 2*dim, 2*dim) = -G.block( 0, 2*dim, 2*dim, 2*dim).conjugate() ;
    }

  return ;

} ;

template void hubbard<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;

template void hubbard<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;

