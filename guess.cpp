#include "common.h"
#include "constants.h"
#include <iostream>
#include <Eigen/Dense>
#include "guess.h"
#include "qtzcntrl.h"
#include "qtzio.h"
#include "solver.h"
#include "util.h"
#include "wfn.h"

void guess_drv( common& com) {
/*
  Given a job to do, generate a guess.  
*/

  return ;

  } ;

void generate_pk( int nbasis, int nele, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa) {
/*
  Generate rho and kappa for an HFB calculation.  Right now there are two things that we can do to generate the wavefunction.
  
    - Thermalize the occupations
    - Mix with a random unitary matrix.
*/

  int f, nbas = 2*nbasis ;
  double temp = 3.0e4 ;
  Eigen::VectorXd A, B ;
  Eigen::MatrixXcd C, D, E ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w ;

  A.resize( nbasis) ;
  B.resize( nbasis) ;
  C.resize( nbas, nbas) ;
  D.resize( nbas, nbas) ;
  E.resize( nbas, nbas) ;
  A.setZero() ;
  B.setZero() ;
  C.setRandom() ;

  w.moc.resize( nbas, nbas) ;
  w.eig.resize( nbas) ;

/*
  Read in a cghf wavefunction
*/
  load_wfn( w) ;

/*
  Generate Thermal occupancies.
*/
  thermal_occv( nbasis, nele, w.eig, temp, A, B) ;

  print_mat( A, " Thermal occ A") ;
  print_mat( B, " Thermal occ B") ;

/*
  Generate a random Unitary matrix
*/
  D = C.adjoint() - C ;
  ahm_exp( D, C, nbas, 0) ;

/*
  Fill U and V according to the Bloch-Messiah decomposition
*/
  D.setZero() ;
  E.setZero() ;

  for( f = 0; f < nbas; f += 2){
    D( f, f+1) = static_cast<cd>(A( f/2)) ;
    D( f+1, f) = static_cast<cd>(-A( f/2)) ;
    E( nbas - f - 1, nbas - f - 1) = static_cast<cd>(B( f/2)) ;
    E( nbas - f - 2, nbas - f - 2) = static_cast<cd>(B( f/2)) ;
    }

  print_mat( D, "V occ") ;
  print_mat( E, "U occ") ;

/*
  use the MO coeffcients and the random mixing to generate U and V
*/
  rho = -C*w.moc*D*D*w.moc.adjoint()*C.adjoint() ;
  kappa = C*w.moc*D*E*w.moc.transpose()*C.transpose() ;

  print_mat( rho, " rho ") ;
  print_mat( kappa, " kappa ") ;

  w.eig.resize( 0) ;
  w.moc.resize( 0, 0) ;
  E.resize( 0, 0) ;
  D.resize( 0, 0) ;
  C.resize( 0, 0) ;
  B.resize( 0) ;
  A.resize( 0) ;

  return ;

  } ;

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
  thermal_occm( nbas, nele, w.eig, temp, A, B) ;

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
  thermal_occm( nbas, nele, w.eig, temp, A, B) ;

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

template< class matrix>
void homo_lumo_mix( matrix& w, const int i) {
/*
  Given a set of MOs, mix the ith and i + 1 orbitals
*/
  int n = w.rows() ;
  typename matrix::Scalar two ;
  two = static_cast<typename matrix::Scalar>( d2) ;
  matrix homo( n, 1), lumo( n, 1), homo_mix( n, 1), lumo_mix( n, 1) ;
/*
  Eigen Matrices start at 0 not 1
*/
  homo = w.col( i-1) ;
  lumo = w.col( i) ;
  homo_mix = (homo + lumo)/two ;
  lumo_mix = (homo - lumo)/two ;
  w.col( i-1) = homo_mix ;
  w.col( i) = lumo_mix ;

  lumo_mix.resize( 0, 0) ;
  homo_mix.resize( 0, 0) ;
  lumo.resize( 0, 0) ;
  homo.resize( 0, 0) ;

  return ;

  } ;

template void homo_lumo_mix( Eigen::MatrixXd&, const int) ;

template void homo_lumo_mix( Eigen::MatrixXcd&, const int) ;

void Uthermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k) {
/*
  Given an initial real Slater Determinant from Hartree-Fock, return
  a thermalized rho and kappa guess mixed with a random Unitary matrix.
*/
  double temp = 3.0e4 ;
  Eigen::MatrixXd A, B ;
  Eigen::MatrixXcd C, D, U ;
  wfn < double, Eigen::Dynamic, Eigen::Dynamic> w ;
  A.resize( nbas, nbas) ;
  B.resize( nbas, nbas) ;
  C.resize( 2*nbas, 2*nbas) ;
  D.resize( 2*nbas, 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;
  A.setZero() ;
  B.setZero() ;
  C.setZero() ;
  D.setZero() ;
  U.setRandom() ;
  w.moc.resize( nbas, nbas) ;
  w.eig.resize( nbas) ;
  load_wfn( w) ;
  thermal_occm( nbas, nele, w.eig, temp, A, B) ;

  D = U.adjoint() - U ;
  D *= z1/z10 ;

  ahm_exp( D, U, 2*nbas, 0) ;

  std::cout << " U.adjoint()*U " << std::endl ;
  std::cout << U.adjoint()*U << std::endl ;

  C.setZero() ;
  D.setZero() ;

  C.block( 0, 0, nbas, nbas) = w.moc ;
  C.block( nbas, nbas, nbas, nbas) = w.moc ;

  D = U*C ;

  print_mat( D, " Mixed Wavefunction ") ;

  C.setZero() ;
  C.block( 0, 0, nbas, nbas).real() = A ;
  C.block( nbas, nbas, nbas, nbas).real() = A ;

  print_mat( C, " Occupations ") ;

  p = D*C*D.adjoint() ;
  print_mat( p, " initial rho ") ;

  C.setZero() ;
  C.block( 0, nbas, nbas, nbas).real() = B ;
  C.block( nbas, 0, nbas, nbas).real() = B ;

  k.setZero() ;
  D.block( 0, nbas, nbas, nbas).setZero() ;
  D.block( nbas, 0, nbas, nbas).setZero() ;

  k = D*C*D.adjoint() ;
  print_mat( k, " initial kappa ") ;

  w.eig.resize( 0) ;
  w.moc.resize( 0, 0) ;
  U.resize( 0, 0) ;
  D.resize( 0, 0) ;
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

void thermal_occm( int& nbas, int& nele, Eigen::Ref<Eigen::VectorXd> eig, double& T, Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) {
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
    B( i, i) = std::sqrt( A( i, i) - A( i, i)*A( i, i)) ;
    }

  return ;

} ;

void thermal_occv( int& nbas, int& nele, Eigen::Ref<Eigen::VectorXd> eig, double& T, Eigen::Ref<Eigen::VectorXd> A, Eigen::Ref<Eigen::VectorXd> B) {
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
      A( i) = d1/( d1 + std::exp( djunk)) ;
      }
    /* nu should be larger than 0 */
    nu = A.sum() - Nele ;

    A.setZero() ;

    for( i = 0; i < nbas; i++){
      djunk = (eig( i) - ef_l)/(kb*T) ;
      A( i) = d1/( d1 + std::exp( djunk)) ;
      }
    /* nl should be smaller than 0 */
    nl = A.sum() - Nele ;

    A.setZero() ;
    
    for( i = 0; i < nbas; i++){
      djunk = (eig( i) - Ef)/(kb*T) ;
      A( i) = d1/( d1 + std::exp( djunk)) ;
      }
    npar = A.sum() - Nele ;

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
    A( i) = d1/( d1 + std::exp( djunk)) ;
    B( i) = std::sqrt( A( i) - A( i)*A( i)) ;
    }

  return ;

} ;
