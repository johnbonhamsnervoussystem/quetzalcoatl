#include "constants.h"
#include<iostream>
#include<string>
#include<vector>
#include<complex>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include "tei.h"
#include "evalm.h"
#include "solver.h"
#include "util.h"

/* This set of routines finds a wavefunction by repeated diagonalization of the fock matrix for 
 * all the flavors of HF.  
 *
 * rrhfdia - real restricted hartree fock
 * crhfdia - complex restricted hartree fock
 * ruhfdia - real unrestricted hartree fock
 * cuhfdia - complex unrestricted hartree fock
 * rghfdia - real generalized hartree fock
 * cghfdia - complex generalized hartree fock
 *
 * */

float rrhfdia( Eigen::Ref<Eigen::MatrixXf> h, Eigen::Ref<Eigen::MatrixXf> s, std::vector<tei>& intarr, int nbasis, int nele, Eigen::Ref<Eigen::MatrixXf> c, Eigen::Ref<Eigen::VectorXf> eig){

  /* Real restricted Hartree-Fock solved by repeated diagonalization. */
  Eigen::MatrixXf f ;
  Eigen::MatrixXf g ;
  Eigen::MatrixXf p ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> f_diag ;
  int iter=0 ;
  int occ ;
  float thresh=1e-8 ;
  float energy ;
  float ene_p=d0 ;
  float e_dif=1e0 ;

  occ = nele/2 ;
  f.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p.resize( nbasis, nbasis) ;

  /* If c has something in it use it as the initial guess. */
  if( c.isZero(0) ) {
    f = h ;
  } else {
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
  } 

  while ( iter < 31 ) {
    iter += 1 ;
    f_diag.compute( f, s) ;
    c = f_diag.eigenvectors().real() ;
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = g ;
    oao ( nbasis, f, s) ;
    f = h + g ;
    g = p*( h + f) ;

    energy = g.trace() ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues() ;
  std::cout << " Density matrix rhf" << std::endl ;
  std::cout << p << std::endl ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  return energy ;

} ;

float crhfdia( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcf> c, Eigen::Ref<Eigen::VectorXf> eig){

  /* Compelx restricted Hartree-Fock solved by repeated diagonalization. */
  typedef std::complex<float> cf ;
  Eigen::MatrixXcf f ;
  Eigen::MatrixXcf g ;
  Eigen::MatrixXcf p ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcf> f_diag ;
  int iter=0 ;
  int occ ;
  float thresh=1e-8 ;
  float energy=d0 ;
  std::complex<float> t_f ;
  float ene_p=d0 ;
  float e_dif=1e0 ;

  occ = nele/2 ;
  f.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p.resize( nbasis, nbasis) ;

  /* If something is saved in c do an inital guess. */
  if( c.isZero(0) ) {
    f = h ;
  } else {
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
  } 

  while ( iter < 31 ) {
    iter += 1 ;
    f_diag.compute( f, s) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
    g = p*(h + f) ;
    t_f = g.trace() ;
    energy = t_f.real() ;
    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues().real() ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  return energy ;

} ;

float ruhfdia( Eigen::Ref<Eigen::MatrixXf> const h, Eigen::Ref<Eigen::MatrixXf> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXf> c_a, Eigen::Ref<Eigen::MatrixXf> c_b, Eigen::Ref<Eigen::VectorXf> eig){

  /* Real unrestricted Hartree-Fock solved by repeated diagonalization. */
  Eigen::MatrixXf f_a ;
  Eigen::MatrixXf f_b ;
  Eigen::MatrixXf g ;
  Eigen::MatrixXf p_a ;
  Eigen::MatrixXf p_b ;
  Eigen::MatrixXf p_t ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> f_diag ;

  int iter=0 ;
  int nbas ;
  float thresh=1e-8 ;
  float energy=d0 ;
  float t_f=d0 ;
  float ene_p=d0 ;
  float e_dif=1e0 ;

  nbas = 2*nbasis ;
  f_a.resize( nbasis, nbasis) ;
  f_b.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p_a.resize( nbasis, nbasis) ;
  p_b.resize( nbasis, nbasis) ;
  p_t.resize( nbasis, nbasis) ;

  /* If c_a and c_b have values use them as the initial guess */
  if ( c_a.isZero(0) ) {
    f_a = h ;
    f_b = h ;
  } else {
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
  } 

  while ( iter < 31 ) {
    iter += 1 ;
    f_diag.compute( f_a, s) ;
    c_a = f_diag.eigenvectors().real() ;
    f_diag.compute( f_b, s) ;
    c_b = f_diag.eigenvectors().real() ;
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    g = p_a*f_a ;
    t_f = g.trace() ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
    g = p_b*f_b ;
    t_f += g.trace() ;
    g = p_t*h ;
    t_f += g.trace() ;
    energy = t_f/2.0 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  f_diag.compute( f_a, s) ;
  eig.head(nbasis) = f_diag.eigenvalues() ;
  f_diag.compute( f_b, s) ;
  eig.tail(nbasis) = f_diag.eigenvalues() ;
  p_t.resize( 0, 0) ;
  p_b.resize( 0, 0) ;
  p_a.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f_b.resize( 0, 0) ;
  f_a.resize( 0, 0) ;

  return energy ;

} ;

float cuhfdia( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXcf> c_a, Eigen::Ref<Eigen::MatrixXcf> c_b, Eigen::Ref<Eigen::VectorXf> eig){

  /* Complex unrestricted Hartree-Fock solved by repeated diagonalization. */
  typedef std::complex<float> cf ;
  Eigen::MatrixXcf f_a ;
  Eigen::MatrixXcf f_b ;
  Eigen::MatrixXcf g ;
  Eigen::MatrixXcf p_a ;
  Eigen::MatrixXcf p_b ;
  Eigen::MatrixXcf p_t ;
  Eigen::MatrixXcf temp ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcf> f_diag ;
  int iter=0 ;
  int nbas ;
  float thresh=1e-8 ;
  float energy=d0 ;
  std::complex<float> t_f ;
  float ene_p=d0 ;
  float e_dif=1e0 ;

  nbas = 2*nbasis ;
  f_a.resize( nbasis, nbasis) ;
  f_b.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p_a.resize( nbasis, nbasis) ;
  p_b.resize( nbasis, nbasis) ;
  p_t.resize( nbasis, nbasis) ;
  temp.resize( nbasis, nbasis) ;

  /* If c_a and c_b are filled use them for an inital guess */
  if ( c_a.isZero(0) ) {
    f_a = h ;
    f_b = h ;
  } else {
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
  } 

  while ( iter < 31 ) {
    iter += 1 ;
    f_diag.compute( f_a, s) ;
    c_a = f_diag.eigenvectors() ;
    f_diag.compute( f_b, s) ;
    c_b = f_diag.eigenvectors() ;
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    g = p_a*f_a ;
    temp = g ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
    g = p_b*f_b ;
    temp = temp + g ;
    g = p_t*h ;
    temp = temp + g ;
    t_f =  temp.trace() ;
    energy = t_f.real()/2.0 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  f_diag.compute( f_a, s) ;
  eig.head(nbasis) = f_diag.eigenvalues().real() ;
  f_diag.compute( f_b, s) ;
  eig.tail(nbasis) = f_diag.eigenvalues().real() ;
  temp.resize( 0, 0) ;
  p_t.resize( 0, 0) ;
  p_b.resize( 0, 0) ;
  p_a.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f_b.resize( 0, 0) ;
  f_a.resize( 0, 0) ;

  return energy ;

} ;

float rghfdia( Eigen::Ref<Eigen::MatrixXf> const h, Eigen::Ref<Eigen::MatrixXf> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXf> c, Eigen::Ref<Eigen::VectorXf> eig){

  /* Real generalized Hartree-Fock solved by repeated diagonalization. */
  Eigen::MatrixXf f ;
  Eigen::MatrixXf g ;
  Eigen::MatrixXf p ;
  Eigen::MatrixXf h_f ;
  Eigen::MatrixXf s_f ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> f_diag ;
  int iter=0 ;
  int nbas ;
  float thresh=1e-8 ;
  float energy=d0 ;
  float ene_p=d0 ;
  float e_dif=1e0 ;

  nbas = nbasis*2 ;
  f.resize( nbas, nbas) ;
  g.resize( nbas, nbas) ;
  p.resize( nbas, nbas) ;
  h_f.resize( nbas, nbas) ;
  s_f.resize( nbas, nbas) ;
  h_f.setZero() ;
  s_f.setZero() ;
  h_f.block( 0, 0, nbasis, nbasis) = h ;
  h_f.block( nbasis, nbasis, nbasis, nbasis) = h ;
  s_f.block( 0, 0, nbasis, nbasis) = s ;
  s_f.block( nbasis, nbasis, nbasis, nbasis) = s ;

  /* If something is stored in c use it for an inital guess */
  if ( c.isZero(0) ) {
    f = h_f ;
  } else {
    c = f_diag.eigenvectors().real() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h_f + g ;
  } 

  while ( iter < 31 ) {
    iter += 1 ;
    f_diag.compute( f, s_f) ;
    c = f_diag.eigenvectors().real() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h_f + g ;
    g = p*( h_f + f) ;
    energy = g.trace()/2.0 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues() ;
  s_f.resize( 0, 0) ;
  h_f.resize( 0, 0) ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  return energy ;

} ;

float cghfdia( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcf> c, Eigen::Ref<Eigen::VectorXf> eig){

  /* Compelx restricted Hartree-Fock solved by repeated diagonalization. */
  typedef std::complex<float> cf ;
  Eigen::MatrixXcf f ;
  Eigen::MatrixXcf g ;
  Eigen::MatrixXcf p ;
  Eigen::MatrixXcf h_f ;
  Eigen::MatrixXcf s_f ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcf> f_diag ;
  int iter=0 ;
  int nbas ;
  float thresh=1e-8 ;
  float energy=d0 ;
  std::complex<float> t_f ;
  float ene_p=d0 ;
  float e_dif=1e0 ;

  nbas = nbasis*2 ;
  f.resize( nbas, nbas) ;
  g.resize( nbas, nbas) ;
  p.resize( nbas, nbas) ;
  h_f.resize( nbas, nbas) ;
  s_f.resize( nbas, nbas) ;
  h_f.setZero() ;
  s_f.setZero() ;
  h_f.block( 0, 0, nbasis, nbasis) = h ;
  h_f.block( nbasis, nbasis, nbasis, nbasis) = h ;
  s_f.block( 0, 0, nbasis, nbasis) = s ;
  s_f.block( nbasis, nbasis, nbasis, nbasis) = s ;

  if ( c.isZero(0) ) {
    f = h_f ;
  } else {
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h_f + g ;
  } 

  while ( iter < 31 ) {
    iter += 1 ;
    f_diag.compute( f, s_f) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h_f + g ;
    g = p*( h_f + f) ;
    t_f = g.trace() ;
    energy = t_f.real()/2.0 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues().real() ;
  s_f.resize( 0, 0) ;
  h_f.resize( 0, 0) ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  return energy ;

} ;

