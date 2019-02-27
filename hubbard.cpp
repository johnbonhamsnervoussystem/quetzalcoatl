#include "constants.h"
#include <Eigen/Core>
#include "evalm.h"
#include "nbodyint.h"
#include <iostream>
#include "qtzio.h"
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
      The routine is written assuming a factor of 1/2 infront of G
      and D.  hubbard does not have those and so we remove them here 
      to maintain the broader implementation.
    */
    G.setZero() ;
    std::cout << " ghfb for hubbard not yet implemented " << std::endl ;
    }

  return ;

} ;

template void hubbard<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;

template void hubbard<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;

