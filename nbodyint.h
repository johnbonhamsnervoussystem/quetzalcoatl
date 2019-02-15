#include <Eigen/Core>

#ifndef NBODYINT_H
#define NBODYINT_H

/*
  At one point I'd like to have this be an eigen matrix object with contraction defined individually for
  each nbody term.  For now however, I will just use it has a container to simplify calling routines.  
*/
template <typename s>
class nbodyint{
  public :
    /*
      Symmetry of the wavefunction - restricted/unrestriced/generalized
    */
    int itype ;
    int dim ;
    Eigen::Matrix<s, Eigen::Dynamic, Eigen::Dynamic> G ;

    nbodyint( int i, int n){
      itype = i ;
      dim = n ;
      G.resize( dim, dim) ;
      } ;

    virtual void contract( Eigen::Ref<Eigen::Matrix<s, Eigen::Dynamic, Eigen::Dynamic>> m){} ;
    Eigen::Matrix<s, Eigen::Dynamic, Eigen::Dynamic> getG( void) ;

    ~nbodyint(){} ;
} ;

#endif
