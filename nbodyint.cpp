#include "constants.h"
#include "nbodyint.h"

template <typename s>
Eigen::Matrix<s, Eigen::Dynamic, Eigen::Dynamic> nbodyint<s>::getG( void){

  return G ;

} ;

template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nbodyint<double>::getG( void) ;
template Eigen::Matrix<cd, Eigen::Dynamic, Eigen::Dynamic> nbodyint<cd>::getG( void) ;

