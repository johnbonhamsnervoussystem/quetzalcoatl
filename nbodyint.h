#include <Eigen/Core>

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
      G.resize( dim, dim) ;
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
