#include "constants.h"
#include <Eigen/Core>
#include "evalm.h"
#include "nbodyint.h"
#include <iostream>
#include "r12.h"

template <typename s>
void r12<s>::contract( Eigen::Ref<Eigen::Matrix<s, Eigen::Dynamic, Eigen::Dynamic>> m){

  if ( itype == 1){
    /*
      rhf
    */
    ctr2er( ints, m, G, dim) ;
  } else if ( itype == 2){
    /*
      uhf
    */
  } else {
    /*
      ghf
    */
    }

  return ;

}

template void r12<double>::contract( Eigen::Ref<Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>> m) ;

template void r12<cd>::contract( Eigen::Ref<Eigen::Matrix< cd, Eigen::Dynamic, Eigen::Dynamic>> m) ;

