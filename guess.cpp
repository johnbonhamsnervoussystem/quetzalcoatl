#include "common.h"
#include "constants.h"
#include <iostream>
#include <Eigen/Dense>
#include "guess.h"
#include "qtzcntrl.h"
#include "qtzio.h"
#include "util.h"
#include "wfn.h"

void guess_drv( common& com) {
/*
  Given a job to do, generate a guess.  
*/

  return ;

}

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k) {
/*
  Given an initial real Slater Determinant from Hartree-Fock, return
  a thermalized rho and kappa guess.
*/
  double temp = 3.0e4 ;
  Eigen::MatrixXd A, B, C ;
  wfn < double, Eigen::Dynamic, Eigen::Dynamic> w ;
  A.resize( nbas, nbas) ;
  B.resize( nbas, nbas) ;
  C.resize( nbas, nbas) ;
  A.setZero() ;
  B.setZero() ;
  C.setZero() ;
  w.moc.resize( nbas, nbas) ;
  w.eig.resize( nbas) ;
  load_wfn( w) ;
//  Ef = (w.eig( nele - 1) + w.eig( nele))/d2 ;
//  Ef = -0.06883 ;
  thermal_occ( nbas, nele, w.eig, temp, A, B) ;
/*
  for( int i=0; i < nbas; i++){
    djunk = (w.eig( i) - Ef)/(kb*3.0e4) ;
    A( i, i) = d1/( d1 + std::exp( djunk)) ;
    B( i, i) = std::sqrt( A( i, i) - A( i, i)*A(i, i)) ;
    }
*/

  C = w.moc*A*w.moc.adjoint() ;
  print_mat( C, " initial rho") ;

  if ( p.rows() == nbas ){
    p.real() = C ;
  } else {
    p.block( 0, 0, nbas, nbas).real() = C ;
    p.block( nbas, nbas, nbas, nbas).real() = C ;
    }

  C = w.moc*B*w.moc.adjoint() ;
  print_mat( C, " initial kappa") ;

  if ( k.rows() == nbas ){
    k.real() = C ;
  } else {
    k.block( 0, nbas, nbas, nbas).real() = C ;
    k.block( nbas, 0, nbas, nbas).real() = -C ;
    }

  w.eig.resize( 0) ;
  w.moc.resize( 0, 0) ;
  C.resize( 0, 0) ;
  B.resize( 0, 0) ;
  A.resize( 0, 0) ;

  return ;

}

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> k) {
/*
  Given an initial real Slater Determinant from Hartree-Fock, return
  a thermalized rho and kappa guess.
*/
  double temp = 3.0e4 ;
  Eigen::MatrixXd A, B, C ;
  wfn < double, Eigen::Dynamic, Eigen::Dynamic> w ;
  A.resize( nbas, nbas) ;
  B.resize( nbas, nbas) ;
  C.resize( nbas, nbas) ;
  A.setZero() ;
  B.setZero() ;
  C.setZero() ;
  w.moc.resize( nbas, nbas) ;
  w.eig.resize( nbas) ;
  load_wfn( w) ;
  thermal_occ( nbas, nele, w.eig, temp, A, B) ;

  C = w.moc*A*w.moc.adjoint() ;
  print_mat( C, " initial rho") ;

  if ( p.rows() == nbas ){
    p = C ;
  } else {
    p.block( 0, 0, nbas, nbas) = C ;
    p.block( nbas, nbas, nbas, nbas) = C ;
    }

  C = w.moc*B*w.moc.adjoint() ;
  print_mat( C, " initial kappa") ;

  if ( k.rows() == nbas ){
    k = C ;
  } else {
    k.block( 0, nbas, nbas, nbas) = C ;
    k.block( nbas, 0, nbas, nbas) = -C ;
    }

  w.eig.resize( 0) ;
  w.moc.resize( 0, 0) ;
  C.resize( 0, 0) ;
  B.resize( 0, 0) ;
  A.resize( 0, 0) ;

  return ;

}

void rand01_guess( Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) {
/*

  A very specific routine that should be generalized.  

  Use a seed and pseudo-random number generator to build identical initial guesses
  with values between -0.5 and 0.5.

*/
  uint32_t seed = 101723 ;
  uint32_t *state ;
  state = &seed ;

  for( int i = 0; i < A.rows(); i++){
    for( int j = 0; j < A.cols(); j++){
      A( i , j) = rand01( state) ;
      }
    }

  for( int i = 0; i < B.rows(); i++){
    for( int j = 0; j < B.cols(); j++){
      B( i , j) = rand01( state) ;
      }
    }

  return ;

}

void thermal_occ( int& nbas, int& nele, Eigen::Ref<Eigen::VectorXd> eig, double& T, Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) {
  /*
    Given a set of eigenvalues and temperature, adjust the fermi energy until the occupation numbers of the density equals nele.
    nbas - number of basis functions
    nele - target number of electrons
    eig - vector of eigenvalues
    T - Temperature to thermalize to
    A - Occupations for density matrix
    B - Occupations for pairing density
  */
  int i, iter = 0 ;
  double nu, nl, npar, djunk, Nele = static_cast<double>(nele) ;
  double ef_u = eig( nele), ef_l = eig( nele - 1), Ef ;

  while ( true){
    if ( iter++ > 500){
      qtzcntrl::shutdown( " Exceeded 500 iterations in thermal_occ ") ;
      return ;
      }

    Ef = ( ef_u + ef_l)/d2 ;

    A.setZero() ;

    for( i = 0; i < nbas; i++){
      djunk = (eig( i) - ef_u)/(kb*T) ;
      A( i, i) = d1/( d1 + std::exp( djunk)) ;
      }
    /* nu should be larger than 0 */
    nu = A.trace() - Nele ;

    A.setZero() ;

    for( i = 0; i < nbas; i++){
      djunk = (eig( i) - ef_l)/(kb*T) ;
      A( i, i) = d1/( d1 + std::exp( djunk)) ;
      }
    /* nl should be smaller than 0 */
    nl = A.trace() - Nele ;

    A.setZero() ;
    
    for( i = 0; i < nbas; i++){
      djunk = (eig( i) - Ef)/(kb*T) ;
      A( i, i) = d1/( d1 + std::exp( djunk)) ;
      }
    npar = A.trace() - Nele ;

    /* Check for convergence */
    if ( std::abs(npar) < 1.0e-5) {
      break ;
    } else if ( nl*nu < d0){
      /* If the boundaries are on opposite sides of particle number continue */
      if ( npar < d0){
        ef_l = Ef ;
        Ef = (ef_l + ef_u)/d2 ;
      } else {
        ef_u = Ef ;
        Ef = (ef_l + ef_u)/d2 ;
        }
    } else {
      /* Expand the boundaries */
      ef_u += 0.1 ;
      ef_l -= 0.1 ;
      }
    }

  std::cout << " Fermi energy found : " << Ef << std::endl ;

  /* We have our fermi energy, populate rho and kappa */
  for( i = 0; i < nbas; i++){
    djunk = ( eig( i) - Ef)/(kb*3.0e4) ;
    A( i, i) = d1/( d1 + std::exp( djunk)) ;
    B( i, i) = std::sqrt( A( i, i) - A( i, i)*A(i, i)) ;
    }

  return ;

} ;
