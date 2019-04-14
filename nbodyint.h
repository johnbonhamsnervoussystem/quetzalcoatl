#include <Eigen/Core>
#include "qtzcntrl.h"

#ifndef NBODYINT_H
#define NBODYINT_H

/*
  At one point I'd like to have this be an eigen matrix object with contraction defined individually for
  each nbody term.  For now however, I will just use it has a container to simplify calling routines.  
*/
template <class matrix>
class nbodyint{
  public :
    /*
      itype = symmetry of the wavefunction - restricted/unrestriced/generalized
      dim   = smallest dimension of the problem
      tr    = symmetric or cannonical orthogonalization
    */
    int itype ;
    int dim ;
    matrix G ;

    nbodyint( int i, int n){
      itype = i ;
      dim = n ;
      switch( itype) {
        case 1 : // rhf
          G.resize( dim, dim) ;
          break ;
        case 2 : // uhf
          G.resize( dim, dim) ;
          break ;
        case 3 : // ghf
          G.resize( 2*dim, 2*dim) ;
          break ;
        case 4 : // rhfb
          G.resize( 2*dim, 2*dim) ;
          break ;
        case 6 : // ghfb
          G.resize( 4*dim, 4*dim) ;
          break ;
        default : // Unimplemented case
          qtzcntrl::shutdown( "Unimplemented case in nbodyint") ;
          break ;
        }
      } ;

    virtual void contract( matrix& m){} ;
    virtual void contract( matrix& m, matrix& n){} ;
    matrix getG( void) ;

    ~nbodyint(){} ;
} ;

#endif
