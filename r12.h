#include <Eigen/Core>
#include "nbodyint.h"
#include <vector>
#include "tei.h"

#ifndef R12_H
#define R12_H


template< typename s>
class r12: public nbodyint<s> {
  private :
    using nbodyint<s>::itype ;
    using nbodyint<s>::dim ;
    using nbodyint<s>::G ;
    std::vector<tei> ints ;
  public :
    r12( std::vector<tei>& intarr, int i, int n) : nbodyint<s>::nbodyint( i, n){
      ints = intarr ;
      } ;
    ~r12(){} ;
    void contract( Eigen::Ref<Eigen::Matrix<s, Eigen::Dynamic, Eigen::Dynamic>> m) ;
} ;

#endif
