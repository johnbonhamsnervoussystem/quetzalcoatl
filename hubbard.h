#include "nbodyint.h"

#ifndef HUBBARD_H
#define HUBBARD_H

template <class matrix>
class hubbard: public nbodyint<matrix> {
  private :
    using nbodyint<matrix>::itype ;
    using nbodyint<matrix>::dim ;
    using nbodyint<matrix>::G ;
    double U ;
    typename matrix::Scalar two ;
  public :
    hubbard( double U, int i, int n) : nbodyint<matrix>::nbodyint( i, n){
      this->U = U ;
      two = static_cast<typename matrix::Scalar>( 2.0e0) ;
      } ;
    ~hubbard(){} ;
    void contract( matrix& m) ;
    void contract( matrix& m, matrix& n) ;
} ;

#endif
