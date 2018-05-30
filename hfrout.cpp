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
  float ene_p=0e0 ;
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

    e_dif = abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues() ;
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
  float energy=0e0 ;
  std::complex<float> t_f ;
  float ene_p=0e0 ;
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
    e_dif = abs(ene_p - energy) ;
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
  float energy=0e0 ;
  float t_f=0e0 ;
  float ene_p=0e0 ;
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

    e_dif = abs(ene_p - energy) ;
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
  float energy=0e0 ;
  std::complex<float> t_f ;
  float ene_p=0e0 ;
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
    std::cout << energy << std::endl ;

    e_dif = abs(ene_p - energy) ;
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
  float energy=0e0 ;
  float ene_p=0e0 ;
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

    e_dif = abs(ene_p - energy) ;
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
  float energy=0e0 ;
  std::complex<float> t_f ;
  float ene_p=0e0 ;
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

    std::cout << energy << std::endl ;
    e_dif = abs(ene_p - energy) ;
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

/* c.row(0)<<cf(0.43387,0.00000),cf(0.27114,-0.00324),cf(-0.09169,-0.57919),cf(-0.38192,-0.20424),cf(-0.26568,-0.16014),cf(-0.26732,-0.21997),cf(-0.24721,-0.24721),cf(-0.04261,0.13549);
 c.row(1)<<cf(0.13605,-0.00368),cf(0.49320,0.00000),cf(0.38220,0.20542),cf(0.09060,0.57876),cf(0.34343,0.29669),cf(0.10034,-0.00656),cf(-0.25158,-0.18272),cf(-0.06120,-0.20488);
 c.row(2)<<cf(0.12976,-0.14974),cf(0.16994,0.11031),cf(-0.04014,0.12340),cf(0.25056,-0.25056),cf(-0.33828,0.11161),cf(0.12925,0.50393),cf(0.36973,0.22312),cf(-0.58322,0.06939);
 c.row(3)<<cf(0.12608,0.14809),cf(0.17318,-0.11176),cf(-0.25037,0.25037),cf(0.04076,-0.12395),cf(0.11855,-0.49750),cf(0.34266,0.13628),cf(0.12906,0.20682),cf(0.68703,0.00000);
 c.row(4)<<cf(0.20081,-0.19883),cf(0.01310,-0.01662),cf(-0.24556,-0.17844),cf(-0.21455,0.06503),cf(0.59865,0.00000),cf(-0.16181,-0.11393),cf(0.67540,0.00000),cf(-0.12687,-0.24329) ;
 c.row(5)<<cf(0.04479,-0.04078),cf(-0.19662,0.19488),cf(-0.21509,0.06471),cf(-0.24500,-0.17868),cf(0.18869,-0.11247),cf(0.59100,0.00000),cf(-0.59327,0.09406),cf(-0.36347,-0.19626) ;
 c.row(6)<<cf(0.35258,-0.18887),cf(-0.09867,0.30340),cf(-0.22339,0.11373),cf(0.68455,0.00000),cf(-0.19553,0.22488),cf(-0.01897,-0.35625),cf(-0.05367,-0.21700),cf(0.23887,0.18807) ;
 c.row(7)<<cf(0.19453,-0.34489),cf(-0.31017,0.09368),cf(0.68403,0.00000),cf(-0.22500,0.11365),cf(-0.08306,-0.18508),cf(-0.13290,0.39656),cf(-0.02846,0.12294),cf(0.25147,0.25147) ; */

/*  c.row(0) << cf(0.573091,0.0), cf(0.0,0.0), cf(0.820822,0.0), cf(0.0,0.0), cf(-0.386664,0.0), cf(0.0,0.0), cf(0.000202978,0.0), cf(0.0,0.0) ;
 c.row(1) << cf(0.117435,0.0), cf(0.0,0.0), cf(0.000179326,0.0), cf(0.0,0.0), cf( 0.681752,0.0), cf(0.0,0.0), cf(0.820436,0.0), cf(0.0,0.0) ;
 c.row(2) << cf(0.117327,0.0), cf(0.0,0.0), cf(0.000166753,0.0), cf(0.0,0.0), cf( 0.680703,0.0), cf(0.0,0.0), cf(-0.821322,0.0), cf(0.0,0.0) ;
 c.row(3) << cf(0.573197,0.0), cf(0.0,0.0), cf(-0.820936,0.0), cf(0.0,0.0), cf(-0.386265,0.0), cf(0.0,0.0), cf(0.000215299,0.0), cf(0.0,0.0) ;
 c.row(4) << cf(0.0,0.0), cf(0.117435,0.0), cf(0.0,0.0), cf(-0.000180394,0.0), cf(0.00,0.0), cf(-0.681748,0.0), cf(0.0,0.0), cf(-0.820439,0.0) ;
 c.row(5) << cf(0.0,0.0), cf(0.573091,0.0), cf(0.0,0.0), cf(-0.820822,0.0), cf(0.00,0.0), cf(0.386665,0.0), cf(0.0,0.0), cf(-0.000201551,0.0) ;
 c.row(6) << cf(0.0,0.0), cf(0.573197,0.0), cf(0.0,0.0), cf(0.820937,0.0), cf(0.00,0.0), cf(0.386264,0.0), cf(0.0,0.0), cf(-0.000214122,0.0) ;
 c.row(7) << cf(0.0,0.0), cf(0.117327,0.0), cf(0.0,0.0), cf(-0.000167762,0.0), cf(0.00,0.0), cf(-0.680706,0.0), cf(0.0,0.0), cf(0.821319,0.0) ; */

/* 
 * These routines check the internal stability of a wavefunction 
 * */
/*  
 *  The Hessian with respect to orbital rotations is 2*occ*virtul square with the 
 *  structure.  i and j refer to the occupied orbitals.  a and b refer to the virtual
 *  orbitals.
 *
 *  [ A   B ]      A = (e_{a} - e_{i})delta_{ia,jb} + (aj||ib)
 *  [ B*  A*]      B = (ab||ij)
 *
 *  (ab||ij) -> (12||12) 
 *
 *  */

/* 
 * Given vectors of the MO coefficients and a list of integrals,
 * contract the MO integrals.
 * */
//float gc_pqrs( std::vector<float> p, std::vector<float> q, std::vector<float> r, std::vector<float> s, intarr, int nb) {
///* Real generalized coulomb contraction of (pq,rs)*?
//    int ntt ;
//    int n2ei ;
//    int i ;
//    int j ;
//    int k ;
//    int l ;
//    float moint=0.0 ;
//    float val ;
//    bool ieqk ;
//    bool jeql ;
//    bool ieqj ;
//    bool keql ;
//
//    ntt = nbasis*(nbasis + 1)/2 ;
//    n2ei = ntt*(ntt+1)/2 ;
//
//    for(int t = 0; t < n2ei; t++ ) {
//      /*  ( 1 1| 2 2)
// *        ( i j| k l) = val */
//      i = intarr[t].r_i() - 1 ;
//      j = intarr[t].r_j() - 1 ;
//      k = intarr[t].r_k() - 1 ;
//      l = intarr[t].r_l() - 1 ;
//      val = intarr[t].r_v() ;
////   
//      ieqj = ( ( i == j) ? true : false) ;
//      keql = ( ( k == l) ? true : false) ;
//      
//      if ( ieqj && keql && i == k ){
//
//        /* ( i i| i i) */
//
//        moint += (p(i)*r(i) + p(i+nb)*r(i+nb))*(q(i)*s(i) + q(i+nb)*s(i+nb))*val ;
//
//      } else if ( not ieqj && keql && j == k ){
//
//        /* ( i j| j j) */
//
//        G( i , j) += p( j, j)*val ;
//        G( j , i) += p( j, j)*val ;
//        G( j , j) += (p( j, i) + p( i, j))*val ;
//
//      } else if ( not ieqj && i == k && j == l ){
//
//        /* ( i j| i j) */
//
//        G( i, j) += (p( j, i) + p( i, j))*val ;
//        G( j, i) += (p( i, j) + p( j, i))*val ;
//
//      } else if ( ieqj && keql && i != k ){
//
//        /* ( i i| k k) */
//
//        G( i, i) += p( k, k)*val ;
//        G( k, k) += p( i, i)*val ;
//
//      } else if ( ieqj && not keql && i == k ){
//
//        /* ( i i| i l) */
//
//        G( i, i) += (p( i, l) + p( l, i))*val ;
//        G( l, i) += p( i, i)*val ;
//        G( i, l) += p( i, i)*val ;
//
//      } else if ( not ieqj && not keql && j == l && i != k){
//
//        /* ( i j| k j) */
//
//        G( i, j) += (p( k, j) + p( j, k))*val ;
//        G( j, i) += (p( k, j) + p( j, k))*val ;
//        G( k, j) += (p( i, j) + p( j, i))*val ;
//        G( j, k) += (p( i, j) + p( j, i))*val ;
//
//      } else if ( not ieqj && keql && i != k && j != k) {
//
//        /* ( i j| k k) */
//
//        G( i, j) += p( k, k)*val ;
//        G( j, i) += p( k, k)*val ;
//        G( k, k) += (p( j, i) + p( i, j))*val ;
//
//      } else if ( ieqj && not keql && i != k && i != l) {
//
//        /* ( i i| k l) */
//
//        G( i, i) += (p( k, l) + p( l, k))*val ;
//        G( k, l) += p( i, i)*val ;
//        G( l, k) += p( i, i)*val ;
//
//      } else if ( not ieqj && not keql && j == k && i != l ) {
//
//        /* ( i j| j l) */
//
//        G( i, j)+= (p( j, l) + p( l, j))*val ;
//        G( j, i)+= (p( j, l) + p( l, j))*val ;
//        G( j, l)+= (p( i, j) + p( j, i))*val ;
//        G( l, j)+= (p( i, j) + p( j, i))*val ;
//
//      } else if ( not ieqj && not keql && i == k && j != l) {
//
//        /* ( i j| i l) */
//
//        G( i, j)+= (p( i, l) + p( l, i))*val ;
//        G( j, i)+= (p( i, l) + p( l, i))*val ;
//        G( i, l)+= (p( i, j) + p( j, i))*val ;
//        G( l, i)+= (p( i, j) + p( j, i))*val ;
//
//      } else {
//
//        /* ( i j| k l) */
//
//        G( i, j)+= (p( k, l) + p( l, k))*val ;
//        G( j, i)+= (p( k, l) + p( l, k))*val ;
//        G( k, l)+= (p( i, j) + p( j, i))*val ;
//        G( l, k)+= (p( i, j) + p( j, i))*val ;
//
//      }
//
//      }
//
//    return 0 ;
//
//  } ;
//}
//
///* The loop goes over occ*virtual 
// * do i -> a
// *      -> b
// *      -> c
// *      ...
// * do j -> a
// *      -> b
// *      -> c
// *      ...
// * etc.
// *
// * Storing this as a two dimentional matix then means
// * 0 = 1st occ, 1st virt
// * nvirt = 2 occ 1st virt
// * 2*nvirt = 3 occ 1st virt
// * etc.
// *
// */
//  for ( int i=0; i < occ; i ++ ){
//    for ( int a=0; a < virt; a++ ){
//      for ( int j=0; j < occ; j ++ ){
//        for ( int b=0; b < virt; b++ ){
//          coul( i,)
//        }
//      }
//    }
//  }


