#include "constants.h"
#include "nbodyint.h"

template < class matrix>
matrix nbodyint<matrix>::getG( void){

  return G ;

} ;

template Eigen::MatrixXd nbodyint<Eigen::MatrixXd>::getG( void) ;
template Eigen::MatrixXcd nbodyint<Eigen::MatrixXcd>::getG( void) ;

/*
template < class matrix>
void nbodyint<matrix>::contract( matrix& m){ return ;} ;

template void nbodyint<Eigen::MatrixXd>::contract( Eigen::MatrixXd&) ;
template void nbodyint<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&) ;

template < class matrix>
void nbodyint<matrix>::contract( matrix& m, matrix& n){ return ;} ;

template void nbodyint<Eigen::MatrixXd>::contract( Eigen::MatrixXd&, Eigen::MatrixXd&) ;
template void nbodyint<Eigen::MatrixXcd>::contract( Eigen::MatrixXcd&, Eigen::MatrixXcd&) ;
*/

