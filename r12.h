#include <Eigen/Core>
#include "nbodyint.h"
#include <vector>
#include "tei.h"

#ifndef R12_H
#define R12_H


template <class matrix>
class r12: public nbodyint<matrix> {
  private :
    using nbodyint<matrix>::itype ;
    using nbodyint<matrix>::dim ;
    using nbodyint<matrix>::G ;
    std::vector<tei>* ints ;
    matrix x ;
    matrix xi ;
  public :
    r12( std::vector<tei>* intarr, matrix& xs, int i, int n) : nbodyint<matrix>::nbodyint( i, n){
      ints = intarr ;
      x.resize( n, n) ;
      xi.resize( n, n) ;
      x = xs ;
      xi = x.inverse() ;
      } ;
    ~r12(){} ;
    void contract( matrix& m) ;
    void contract( matrix& m, matrix& n) ;
} ;

#endif
