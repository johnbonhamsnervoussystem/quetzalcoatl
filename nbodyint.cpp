#include "constants.h"
#include "nbodyint.h"

template < class matrix>
matrix nbodyint<matrix>::getG( void){

  return G ;

} ;

template Eigen::MatrixXd nbodyint<Eigen::MatrixXd>::getG( void) ;
template Eigen::MatrixXcd nbodyint<Eigen::MatrixXcd>::getG( void) ;

