#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include <complex>
#include <math.h>
#include "wigner.h"
#include "common.h"
#include "util.h"

double d_cof ( int j, int m, int k, int n) {
/*
 * Calculate the denominator for the prefactor since the numerator has no 'n' 
 * dependence.
 *
 * This assumes the integers given will not produce negative factorials.
 *
 * 1/[(j +k -n)!(j -m -n)!n!( n +m -k)!]
 *
 */
  double val=1.0 ;
  double cof_d=1.0 ;

// Do the denominator

  cof_d = cof_d*fact(j + k - n) ;
  cof_d = cof_d*fact(j - m - n) ;
  cof_d = cof_d*fact( n) ;
  cof_d = cof_d*fact(n + m - k) ;

  val = val/cof_d ;

  return val ;

} ;

double d_cof ( double j, double m, double k, int n) {
/*
 * Calculate the denominator for the prefactor since the numerator has no 'n' 
 * dependence.
 *
 * This assumes the integers given will not produce negative factorials.
 *
 * 1/[(j +k -n)!(j -m -n)!n!( n +m -k)!]
 *
 */
  double val=1.0 ;
  double cof_d=1.0 ;

// Do the denominator

  cof_d = cof_d*fact( static_cast<int>( j + k) - n) ;
  cof_d = cof_d*fact( static_cast<int>( j - m) - n) ;
  cof_d = cof_d*fact( n) ;
  cof_d = cof_d*fact( n + static_cast<int>( m - k)) ;

  val = val/cof_d ;

  return val ;

} ;

double small_wd ( int j, int m, int k, double beta) {
/*
 * Calculate 
 * d^{j}_{mk}(Theta) = sum_{n_min}^[n_max}(-1)^{n}W_{n}^{jmk}
 * for integer values of j, m, k
 */ 
  int cosexp ;
  int sinexp ;
  int n_min ;
  int n_max ;
  double cos_b ;
  double sin_b ;
  double cof_n = 1.0 ;
  double d_jmk = d0 ;
  double w_n ;
 
  n_min = std::max( 0, k - m) ;
  n_max = std::min( j - m, j + k) ;

/* 
 * The numerator of the prefactor has no dependence on n so 
 * we calculate it once.
 * [(j +m)!(j -m)!(j +k)!(j -k)!]^{1/2} 
 * */

  cof_n = cof_n*fact(j + m) ;
  cof_n = cof_n*fact(j - m) ;
  cof_n = cof_n*fact(j + k) ;
  cof_n = cof_n*fact(j - k) ;
  cof_n = sqrt( cof_n) ;

  for ( int n = n_min ; n <= n_max ; n++) {
    cosexp = 2*j + k - m - 2*n ;
    sinexp = m - k + 2*n ;
    sin_b = sin( beta/2.0) ;
    cos_b = cos( beta/2.0) ;
    w_n = d_cof( j, m, k, n)*pow( cos_b, cosexp)*pow( sin_b, sinexp) ;
    d_jmk += w_n*pow( -1.0, n) ;
  }

  d_jmk = d_jmk*cof_n ;

  return d_jmk ;

} ;

double small_wd ( double j, double m, double k, double beta) {
/*
 * Calculate 
 * d^{j}_{mk}(Theta) = sum_{n_min}^[n_max}(-1)^{n}W_{n}^{jmk}
 * for integer values of j, m, k
 */ 
  int cosexp ;
  int sinexp ;
  int n_min ;
  int n_max ;
  double cos_b ;
  double sin_b ;
  double cof_n = 1.0 ;
  double d_jmk = d0 ;
  double w_n ;
 
  n_min = std::max( 0, static_cast<int>(k - m)) ;
  n_max = std::min( static_cast<int>(j - m), static_cast<int>(j + k)) ;

/* 
 * Overloaded version of small_wd to handle half integer values.  
 * I'd be very surprised to find out there is not a better method of
 * converting these doubles to integers but I'll chalk it up to
 * my c++ skills for now.
 * */

  cof_n = cof_n*fact( static_cast<int>(j + m)) ;
  cof_n = cof_n*fact( static_cast<int>(j - m)) ;
  cof_n = cof_n*fact( static_cast<int>(j + k)) ;
  cof_n = cof_n*fact( static_cast<int>(j - k)) ;
  cof_n = sqrt( cof_n) ;

  for ( int n = n_min ; n <= n_max ; n++) {
    cosexp = static_cast<int>(2.0*j) + static_cast<int>(k - m) - 2*n ;
    sinexp = static_cast<int>(m - k) + 2*n ;
    sin_b = sin( beta/2.0) ;
    cos_b = cos( beta/2.0) ;
    w_n = d_cof( j, m, k, n)*pow( cos_b, cosexp)*pow( sin_b, sinexp) ;
    d_jmk += w_n*pow( -1.0, n) ;
  }

  d_jmk = d_jmk*cof_n ;

  return d_jmk ;

} ;

std::complex<double> wigner_D ( int j, int m, int k, double alpha, double beta, double gamma) {
/*
  D^(j)_{m',m}( alpha, beta, gamma) = Exp[-im' alpha]*d^(j)_{m',m}( beta)*Exp[-im gamma]
 */
  double d_jmk ;
  cd w_D ;
  cd exp_m ;
  cd exp_k ;
  cd m_arg ;
  cd k_arg ;
  const cd di(0.0, 1.0) ;
  
  d_jmk = small_wd ( j, m, k, beta) ;
  m_arg = -di*cd( m, d0)*cd( alpha, d0) ;
  k_arg = -di*cd( k, d0)*cd( gamma, d0) ;
  exp_m = std::exp( m_arg) ;
  exp_k = std::exp( k_arg) ;
  w_D = exp_m*cd( d_jmk, d0)*exp_k ;

  return w_D ;

} ;

std::complex<double> wigner_D ( double j, double m, double k, double alpha, double beta, double gamma) {

/*
 * D^(j)_{m,k}( alpha, beta, gamma) = Exp[-im alpha]*d^(j)_{m ,k}( beta)*Exp[-ik gamma]
 */

  double d_jmk ;
  cd w_D ;
  cd exp_m ;
  cd exp_k ;
  cd m_arg ;
  cd k_arg ;
  const cd di(0.0, 1.0) ;
  
  d_jmk = small_wd ( j, m, k, beta) ;
  m_arg = di*cd( m, d0)*cd( alpha, d0) ;
  k_arg = di*cd( k, d0)*cd( gamma, d0) ;
  exp_m = std::exp( m_arg) ;
  exp_k = std::exp( k_arg) ;
  w_D = exp_m*cd( d_jmk, d0)*exp_k ;

  return w_D ;

} ;

//void R_s ( common& c, hfwfn& a, hfwfn& b, double alpha, double beta, double gamma) {
//
///* This routine does a spin rotation on a determinant using Euler angles
// * alpha, beta and gamma 
// *
// * Exp[-alpha S_{z}]*Exp[-beta S_{y}]*Exp[-gamma S_{z}]|det>
// *
// * a - input hartree-fock wavefunction
// * b - output hartree-fock wavefunction
// * alpha - angle of the S_{z} rotation
// * beta - angle of the S_{x} rotation
// * gamma - angle of the S_{z} rotation
// *
// * */
//
//  int nbas ;
//  cd Cga ;
//  cd zega ;
//  cd b_cos ;
//  cd b_sin ;
//  const cd di(0.0, 1.0) ;
//  Eigen::MatrixXcd moa ;
//  Eigen::MatrixXcd mob ;
//
//  nbas = c.nbas() ;
//
//  moa.resize( 2*nbas, 2*nbas) ;
//  mob.resize( 2*nbas, 2*nbas) ;
//  mob.setZero() ;
//
//  a.get_mos( moa) ;
//  /* Gamma rotation of the alpha block */
//  zega = -cd( d0, gamma/2.0) ;
//  Cga = std::exp( zega) ;
//  moa.block( 0, 0, nbas, 2*nbas) = Cga*moa.block( 0, 0, nbas, 2*nbas) ;
// 
//  /* Gamma rotation of the beta block */
//  zega = cd( d0, gamma/2.0) ;
//  Cga = std::exp( zega) ;
//  moa.block( nbas, 0, nbas, 2*nbas) = Cga*moa.block( nbas, 0, nbas, 2*nbas) ;
//
//  /* beta rotation mixing */
//  b_cos = cd( cos(beta/2.0), d0) ;
//  b_sin = cd( sin(beta/2.0), d0) ;
//  // Alpha Block
//  mob.block( 0, 0, nbas, 2*nbas) = b_cos*moa.block( 0, 0, nbas, 2*nbas) - b_sin*moa.block( nbas, 0, nbas, 2*nbas) ;
// // Beta Block
//  mob.block( nbas, 0, nbas, 2*nbas) = b_cos*moa.block( nbas, 0, nbas, 2*nbas) + b_sin*moa.block( 0, 0, nbas, 2*nbas) ;
//
//  /* Alpha rotation of the alpha block */
//  zega = -cd( d0, alpha/2.0) ;
//  Cga = std::exp( zega) ;
//  mob.block( 0, 0, nbas, 2*nbas) = Cga*mob.block( 0, 0, nbas, 2*nbas) ;
//
//  /* Alpha rotation of the beta block */
//  zega = cd( d0, alpha/2.0) ;
//  Cga = std::exp( zega) ;
//  mob.block( nbas, 0, nbas, 2*nbas) = Cga*mob.block( nbas, 0, nbas, 2*nbas) ;
//
//  b.set_mos( mob) ;
//  mob.resize( 0, 0) ;
//  moa.resize( 0, 0) ;
//
//  return ;
//
//} ;

//void R_s ( int nbas, hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mob, double alpha, double beta, double gamma) {
//
///* This routine does a spin rotation on a determinant using Euler angles
// * alpha, beta and gamma 
// *
// * Exp[-alpha S_{z}]*Exp[-beta S_{y}]*Exp[-gamma S_{z}]|det>
// *
// * a - input hartree-fock wavefunction
// * b - output hartree-fock wavefunction
// * alpha - angle of the S_{z} rotation
// * beta - angle of the S_{x} rotation
// * gamma - angle of the S_{z} rotation
// *
// * */
//
//  cd Cga ;
//  cd zega ;
//  cd b_cos ;
//  cd b_sin ;
//  const cd di(0.0, 1.0) ;
//  Eigen::MatrixXcd moa ;
//
//  moa.resize( 2*nbas, 2*nbas) ;
//  mob.setZero() ;
//
//  a.get_mos( moa) ;
//  /* Gamma rotation of the alpha block */
//  zega = -cd( d0, gamma/2.0) ;
//  Cga = std::exp( zega) ;
//  moa.block( 0, 0, nbas, 2*nbas) = Cga*moa.block( 0, 0, nbas, 2*nbas) ;
// 
//  /* Gamma rotation of the beta block */
//  zega = cd( d0, gamma/2.0) ;
//  Cga = std::exp( zega) ;
//  moa.block( nbas, 0, nbas, 2*nbas) = Cga*moa.block( nbas, 0, nbas, 2*nbas) ;
//
//  /* beta rotation mixing */
//  b_cos = cd( cos(beta/2.0), d0) ;
//  b_sin = cd( sin(beta/2.0), d0) ;
//  // Alpha Block
//  mob.block( 0, 0, nbas, 2*nbas) = b_cos*moa.block( 0, 0, nbas, 2*nbas) - b_sin*moa.block( nbas, 0, nbas, 2*nbas) ;
// // Beta Block
//  mob.block( nbas, 0, nbas, 2*nbas) = b_cos*moa.block( nbas, 0, nbas, 2*nbas) + b_sin*moa.block( 0, 0, nbas, 2*nbas) ;
//
//  /* Alpha rotation of the alpha block */
//  zega = -cd( d0, alpha/2.0) ;
//  Cga = std::exp( zega) ;
//  mob.block( 0, 0, nbas, 2*nbas) = Cga*mob.block( 0, 0, nbas, 2*nbas) ;
//
//  /* Alpha rotation of the beta block */
//  zega = cd( d0, alpha/2.0) ;
//  Cga = std::exp( zega) ;
//  mob.block( nbas, 0, nbas, 2*nbas) = Cga*mob.block( nbas, 0, nbas, 2*nbas) ;
//
//  moa.resize( 0, 0) ;
//
//  return ;
//
//} ;

void R_s ( int nbas, Eigen::Ref<Eigen::MatrixXcd> moa, Eigen::Ref<Eigen::MatrixXcd> mob, double alpha, double beta, double gamma) {

/* This routine does a spin rotation on a determinant using Euler angles
 * alpha, beta and gamma 
 *
 * Exp[-alpha S_{z}]*Exp[-beta S_{y}]*Exp[-gamma S_{z}]|det>
 *
 * a - input hartree-fock wavefunction
 * b - output hartree-fock wavefunction
 * alpha - angle of the S_{z} rotation
 * beta - angle of the S_{x} rotation
 * gamma - angle of the S_{z} rotation
 *
 * */

  cd Cga ;
  cd zega ;
  cd b_cos ;
  cd b_sin ;
  const cd di(0.0, 1.0) ;

  mob.setZero() ;

  /* Gamma rotation of the alpha block */
  zega = -cd( d0, gamma/2.0) ;
  Cga = std::exp( zega) ;
  moa.block( 0, 0, nbas, 2*nbas) = Cga*moa.block( 0, 0, nbas, 2*nbas) ;
 
  /* Gamma rotation of the beta block */
  zega = cd( d0, gamma/2.0) ;
  Cga = std::exp( zega) ;
  moa.block( nbas, 0, nbas, 2*nbas) = Cga*moa.block( nbas, 0, nbas, 2*nbas) ;

  /* beta rotation mixing */
  b_cos = cd( cos(beta/2.0), d0) ;
  b_sin = cd( sin(beta/2.0), d0) ;
  // Alpha Block
  mob.block( 0, 0, nbas, 2*nbas) = b_cos*moa.block( 0, 0, nbas, 2*nbas) - b_sin*moa.block( nbas, 0, nbas, 2*nbas) ;
 // Beta Block
  mob.block( nbas, 0, nbas, 2*nbas) = b_cos*moa.block( nbas, 0, nbas, 2*nbas) + b_sin*moa.block( 0, 0, nbas, 2*nbas) ;

  /* Alpha rotation of the alpha block */
  zega = -cd( d0, alpha/2.0) ;
  Cga = std::exp( zega) ;
  mob.block( 0, 0, nbas, 2*nbas) = Cga*mob.block( 0, 0, nbas, 2*nbas) ;

  /* Alpha rotation of the beta block */
  zega = cd( d0, alpha/2.0) ;
  Cga = std::exp( zega) ;
  mob.block( nbas, 0, nbas, 2*nbas) = Cga*mob.block( nbas, 0, nbas, 2*nbas) ;

  return ;

} ;

