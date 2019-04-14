#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "binio.h"
#include "common.h"
#include "solver.h"
#include "tei.h"
#include "util.h"
#include "evalm.h"
#include "qtzio.h"

/*
 * Notes for Solver:
 *
 *   - KT 06/2018 Added ahm_exp here because it makes use of the 
 *     Eigenvalues package. 
 *   - Find a real convergence criteria for ahm_exp based on repreated
 *     mutliplication
 *
 * */

void symort( Eigen::Ref<Eigen::MatrixXd> const s, Eigen::Ref<Eigen::MatrixXcd> T) {

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s_diag( s) ;
  /*
    Square matrices only please
  */
  Eigen::MatrixXd::Index dim = s.rows() ;
  Eigen::MatrixXcd s_U ( dim, dim) ;
  Eigen::VectorXd s_eig (dim) ;
  
  s_U = s_diag.eigenvectors() ;
  s_eig = s_diag.eigenvalues() ;
  
  s_eig = d1/sqrt( s_eig.array()) ;
  T = s_U*s_eig.asDiagonal()*s_U.adjoint() ;

  s_eig.resize( 0) ;
  s_U.resize( 0, 0) ;

  return ;

} ;

void canort( Eigen::Ref<Eigen::MatrixXd> const s, Eigen::Ref<Eigen::MatrixXcd> xs) {
  /*  
 *    Find an orthogonal transformation using canonical orthogonalization for a real matrix
 *    */

  int indx=-1 ;
  double thresh = 1e-7 ;
  Eigen::MatrixXd::Index dim = s.rows() ;
  std::complex<double> eig ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s_diag( s) ;
  Eigen::MatrixXcd s_U( dim, dim) ;
  Eigen::VectorXcd s_eig( dim) ;

  /* Allocate some space to build our transformation matrix. */

  s_eig = s_diag.eigenvalues() ;

  s_U = s_diag.eigenvectors() ;

  xs = s_U.householderQr().householderQ() ;

  for (int i = 0; i < dim; i++ ) {
    eig = s_eig(i) ;
    if ( std::real(eig) <= thresh && std::imag(eig) <= thresh ) {
      std::cout << "Linear dependency found in canort.  Reducing basis." << std::endl ;
      std::cout << s_eig(i) << " is less than the threshhold" << std::endl ;
      std::cout << "Reducing basis" << std::endl ;
    } else {
      indx ++ ;
    }
    xs.col(indx) = xs.col(i)/std::sqrt(eig) ;
  }
 
  if ( indx != dim-1 ){
    xs.resize(dim,indx) ;
  }

  dim = indx + 1 ;
 
  /* Deallocate  */

  s_eig.resize( 0) ;
  s_U.resize( 0, 0) ;

  return ;

} ;

void canort( Eigen::Ref<Eigen::MatrixXcd> const s, Eigen::Ref<Eigen::MatrixXcd> xs) {
  /*  
 *    Find an orthogonal transformation using canonical orthogonalization for a complex matrix
 *    */

  int indx=-1 ;
  double thresh = 1.0e-7 ;
  Eigen::MatrixXcd::Index dim = s.rows() ;
  std::complex<double> eig ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> s_diag( s) ;
  Eigen::MatrixXcd s_U ( dim, dim) ;
  Eigen::VectorXcd s_eig( dim) ;

  /* Allocate some space to build our transformation matrix. */

  s_eig = s_diag.eigenvalues() ;

  s_U = s_diag.eigenvectors() ;

  xs = s_U.householderQr().householderQ() ;

  for (int i = 0; i < dim; i++ ) {
    eig = s_eig(i) ;
    if ( std::real(std::abs(eig)) <= thresh ) {
      std::cout << "Linear dependency found in canort.  Reducing basis." << std::endl ;
      std::cout << s_eig(i) << " is less than the threshhold" << std::endl ;
      std::cout << "Reducing basis" << std::endl ;
    } else {
      indx ++ ;
    }
    xs.col(indx) = xs.col(i)/std::sqrt(eig) ;
  }
 
//  if ( indx != dim-1 ){
//    xs.resize(4,indx) ;
//  }

  dim = indx + 1 ;
 
  /* Deallocate  */

  s_eig.resize( 0) ;
  s_U.resize( 0, 0) ;

  return ;

} ;

void ahm_exp( Eigen::Ref<Eigen::MatrixXcd> x, Eigen::Ref<Eigen::MatrixXcd> u, int dim, int opt) {

/*
  Anti-Hermetian Matrix EXPonential

  Given an anti-hermetian matrix in an exponential, build a unitary
  matrix. 

  This is done two ways. The first is by diagonalizing the 
  anti-hermetian matrix

    X = iVdVt     VtV = 1

  U = exp(X) = exp(iVdVt) = V*exp(id)*Vt

  The other method is by the definition of a matrix exponential

  U = exp(X) = 1 +  sum_{n=1}^{infty} X^{n}/n!

  where 1 is the unit matrix.
*/

  int n = 0 ;
  double denom=1.0 ;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> xs ;
  Eigen::MatrixXcd d ;
  Eigen::MatrixXcd tmp ;

  tmp.resize( dim, dim) ;
  d.resize( dim, dim) ;
  d.setIdentity() ;

  if ( opt == 0 ) {
    /* Default behavior.  Get the eigenvalues and vectors of x */
    xs.compute( x) ;

    /* Build d */
    for ( int i = 0; i < dim; i++ ){
      d( i, i) = std::exp(xs.eigenvalues()[i]) ;
    }

    /* Build U by Vtexp(idelta)V */
    tmp = d*xs.eigenvectors() ;
    u = xs.eigenvectors().adjoint()*tmp ;

  } else if ( opt == 1 ) {
    /* Find U by repeated multiplication */
      u.setZero() ;

      /* What should the convergence criteria be? */
      while( n < 15 ) {
        u += d ;
        n ++ ;
        tmp = d*x ;
        denom = 1.0/static_cast<double>(n) ;
        d = denom*tmp ;
      }

  } else {

    std::cout << "Unrecgonized opt in ahm_exp" << std::endl ;
    std::cout << "Returning Unit matrix" << std::endl ;
    u.setIdentity() ;

  }

  return ;

} ;

