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
      Symmetry of the wavefunction - restricted/unrestriced/generalized
    */
    int itype ;
    int dim ;
    matrix G ;

    nbodyint( int i, int n){
      itype = i ;
      dim = n ;
      switch( itype) {
        case 1 : // rrhf
          G.resize( dim, dim) ;
          break ;
        case 2 : // ruhf
          G.resize( dim, dim) ;
          break ;
        case 3 : // rghf
          G.resize( 2*dim, 2*dim) ;
          break ;
        default : // Unimplemented case
          qtzcntrl::shutdown( "Unimplemented case in nbodyint") ;
          break ;
        }
      } ;

/*
    void contract( matrix& m) ;
    void contract( matrix& m, matrix& n) ;
*/

    virtual void contract( matrix& m){} ;
    virtual void contract( matrix& m, matrix& n){} ;

/*
    virtual void contract( Eigen::MatrixXd& m){} ;
    virtual void contract( Eigen::MatrixXcd& m){} ;
    virtual void contract( Eigen::MatrixXd& m, Eigen::MatrixXd& n){} ;
    virtual void contract( Eigen::MatrixXcd& m, Eigen::MatrixXcd& n){} ;
*/
    matrix getG( void) ;

    ~nbodyint(){} ;
} ;

#endif
