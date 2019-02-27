#include <Eigen/Core>
#include <iostream>
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
    typename matrix::Scalar two ;
  public :
    r12( std::vector<tei>* intarr, matrix& xs, int i, int n) : nbodyint<matrix>::nbodyint( i, n){
      ints = intarr ;
      two = static_cast<typename matrix::Scalar>( 2.0e0) ;
      if ( itype == 6){
        x.resize( 2*n, 2*n) ;
        xi.resize( 2*n, 2*n) ;
      } else {
        x.resize( n, n) ;
        xi.resize( n, n) ;
        }
      x = xs ;
      xi = x.inverse() ;
      } ;
    ~r12(){} ;
    void contract( matrix& m) ;
    void contract( matrix& m, matrix& n) ;
} ;

#endif
