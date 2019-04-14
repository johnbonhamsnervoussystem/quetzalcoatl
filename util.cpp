#include <algorithm>
#include "basis.h"
#include <cmath>
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_hyperg.h>
#include "integr.h"
#include <iostream>
#include <Eigen/LU>
#include <Eigen/QR>
#include "qtzio.h"
#include "solver.h"
#include "tei.h"
#include <time.h>
#include "time_dbg.h"
#include "util.h"
#include <vector>
#include "wigner.h"
#include "wfn.h"

/* Utilities that don't belong elsewhere. */

void bonnet_r ( int n, double x, double& p2) {
/* 
  Use Bonnet's recursion formula to evaluate the nth Legendre Polynomial
  at the point x.
 
  P_{m+1}(x) = ((2n + 1)xP_{m}(x) - nP_{n-1}(x))/(m + 1)
*/
  double p0 ;
  double p1 = 1.0 ;

  // Handle edge cases and
  if ( n < 0 ){
   std::cout << "Warning: Negative polynomial in bonnet_r" << std::endl ;
   return ;
  }

  if ( n == 0 ){ p2 = 1.0; return ; }
  if ( n == 1 ){ p2 = x ; return ; }

  p2 = x ;
 
  for ( int i = 1; i < n ; i++){
    p0 = p1 ;
    p1 = p2 ;
    p2 = (static_cast<double>(2*i + 1)*x*p1 - i*p0)/static_cast<double>(i + 1) ;
  }
 
  return ;

  }
  
void bonnet_r( int n, double x, double& p2, double& dp) {
/* 
  Use Bonnet's recursion formula to evaluate the nth Legendre Polynomial
  at the point x.  Additionally, return the derivative of P_{n}
 
  P_{m+1}(x) = ((2n + 1)xP_{m}(x) - nP_{n-1}(x))/(m + 1)
*/ 
  double p0 ;
  double p1 = 1.0 ;

  // Handle edge cases and
  if ( n < 0 ){
   std::cout << "Warning: Negative polynomial in bonnet_r" << std::endl ;
   return ;
  }

  if ( n == 0 ){ p2 = 1.0; dp=d0; return ; }
  if ( n == 1 ){ p2 = x; dp=1.0; return ; }

  p2 = x ;
 
  for ( int i = 1; i < n ; i++){
    p0 = p1 ;
    p1 = p2 ;
    p2 = (static_cast<double>(2*i + 1)*x*p1 - i*p0)/static_cast<double>(i + 1) ;
  }

  dp = static_cast<double>(n)*(x*p2 - p1)/(x*x - 1.0) ;

  return ;

  }

void charge_magnetic_decomposition( const Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> mx, Eigen::Ref<Eigen::MatrixXcd> my, Eigen::Ref<Eigen::MatrixXcd> mz){

/*

  Given a density matrix of the form

    | paa  pab|
    | pba  pbb|

  c = ( paa + pbb)/2
  mx = ( pba + pab)/2
  my = i( pab - pba)/2
  mz = ( paa - pbb)/2

  Decompose it into charge and magnetization density matrices

*/

  int n = c.cols() ;
  
  c = (p.block( 0, 0, n, n) + p.block( n, n, n, n))/z2 ;
  mx = (p.block( 0, n, n, n) + p.block( n, 0, n, n))/z2 ;
  my = zi*(p.block( 0, n, n, n) - p.block( n, 0, n, n))/z2 ;
  mz = (p.block( 0, 0, n, n) - p.block( n, n, n, n))/z2 ;

  return ;

  } ;


  double fact ( int n) {
    /* Return the factorial of n */
    double f = d1 ;
  
    if ( n <= 1 ) { return f ;}
  
    while ( n > 1 ) {
      f = f*static_cast<double>(n) ;
      n+= -1 ;
    }
  
    return f ;
  
  } ;

  double factfact ( int n) {
    /* Return the double factorial of n */
    double f=1.0 ;
  
    if ( n < 1 ) { return 1.0 ;}
  
    while ( n > 1 ) {
      f = f*double(n) ;
      n+= -2 ;
    }
  
    return f ;
  
  } ;

double fboys( int k, double t) {
    double f, m, n ;
/*
Let's not mess around.  Just compute this
*/

    m = static_cast<double>(k) + d1/d2 ;
    n = static_cast<double>(k) + d3/d2 ;
    f = gsl_sf_hyperg_1F1( m, n, -t)/( d2*static_cast<double>(k) + d1) ;
    
    return f ;

  } ;

void fc_hamiltonian ( Eigen::Ref<Eigen::MatrixXcd> H_fc, basis_set b, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::MatrixXd> ns) {

  int n = b.b.size(), natm = c.rows() ;
  int i, j, k, l, m ;
  int q ;
  double t, u, x, y, z, r ;
  atm_basis ai, ak ;
  sto zj, zl ;
  Eigen::MatrixXcd I ( n, n) ;
  I.setConstant( z1) ;

  
  /*
    Loop over the basis set 
  */
  /*
    Atoms first
  */
  for ( i = 0; i < n; i++){
    ai = b.b[i] ;
    /*
      stos on that atom
    */
    for ( j = 0; j < ai.nshl; j++){
      zj = ai.s[j] ;
      for ( k = 0; k < n; k++){
        ak = b.b[k] ;
        for ( l = 0; l < ak.nshl; l++){
          zl = ak.s[l] ;
          /*
            Loop over the nuclear coordinates
          */
          for ( m = 0; m < natm; m++){
            /*
              Given the coordinates of basis functions r and s and the nucleus,
              generate the value of phi(R)_{r}^{*}phi(R)_{s}^{*}
            */
            x = c( m, 0) - ai.c[0] ;
            y = c( m, 1) - ai.c[1] ;
            z = c( m, 2) - ai.c[2] ;
            r = x*x + y*y + z*z ;
            /*
              Iterate over the primitives of the sto
            */
            t = 0.0 ;
            for ( q = 0; q < zj.nprm; q++){
              t += zj.g[q].c*std::pow( x, zj.l[0])*std::pow( y, zj.l[1])*std::pow( z, zj.l[2])*std::exp( -zj.g[q].x*r)/zj.norm ;
              }

            x = c( m, 0) - ak.c[0] ;
            y = c( m, 1) - ak.c[1] ;
            z = c( m, 2) - ak.c[2] ;
            r = x*x + y*y + z*z ;

            u = 0.0 ;
            for ( q = 0; q < zl.nprm; q++){
              u += zl.g[q].c*std::pow( x, zl.l[0])*std::pow( y, zl.l[1])*std::pow( z, zl.l[2])*std::exp( -zl.g[q].x*r)/zl.norm ;
              }
            /*
              Now that we have the magnitude, add each component to the FC Hamiltonian Perturbation term
            */
            H_fc.block( 0, 0, n, n) += I*static_cast<cd>(u*t)*static_cast<cd>(ns( m, 2))/z2 ;
            H_fc.block( n, n, n, n) -= I*static_cast<cd>(u*t)*static_cast<cd>(ns( m, 2))/z2 ;
            H_fc.block( 0, n, n, n) += I*static_cast<cd>(u*t)*(static_cast<cd>(ns( m, 0)) + zi*static_cast<cd>(ns( m, 1)))/z2 ;
            H_fc.block( n, 0, n, n) += I*static_cast<cd>(u*t)*(static_cast<cd>(ns( m, 0)) - zi*static_cast<cd>(ns( m, 1)))/z2 ;
            }
          }
        }
      }
    }

  return ;

  } ;

void transform( int o, Eigen::Ref<Eigen::MatrixXd> X, Eigen::Ref<Eigen::MatrixXd> M) {
/*
  Given a tranformation matrix, convert the passed array to or from the 
  basis transformation passed.  For now we will assume the passed matrices
  are square.  

  options :
    o = 0 - Perform X^{t}*M*X
    o = 1 - Perform X*M*X^{t}
    o = 2 - Perform X*M*X
    o = 3 - Perform X^{t}*M*X^{t}
    o = 4 - Perform X^{T}*M*X
*/

  Eigen::MatrixXd::Index xcol = X.cols(), mcol = M.cols() ;
  Eigen::MatrixXd T ( xcol, xcol) ;
  
  if ( xcol == mcol){
    if ( o == 0){
      T = X.adjoint()*M ;
      M = T*X ;
    } else if ( o == 1){
      T = X*M ;
      M = T*X.adjoint() ;
    } else if ( o == 2){
      T = X*M ;
      M = T*X ;
    } else if ( o == 3){
      T = X.adjoint()*M ;
      M = T*X.adjoint() ;
    } else if ( o == 4){
      T = X.transpose()*M ;
      M = T*X ;
      }
  } else {
    if ( o == 0){
      T = X.adjoint()*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X ;
      T = X.adjoint()*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X ;
      T = X.adjoint()*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X ;
      T = X.adjoint()*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X ;
    } else if ( o == 1){
      T = X*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X.adjoint() ;
      T = X*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X.adjoint() ;
      T = X*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X.adjoint() ;
      T = X*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X.adjoint() ;
    } else if ( o == 2){
      T = X*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X ;
      T = X*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X ;
      T = X*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X ;
      T = X*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X ;
    } else if ( o == 3){
      T = X.adjoint()*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X.adjoint() ;
      T = X.adjoint()*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X.adjoint() ;
      T = X.adjoint()*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X.adjoint() ;
      T = X.adjoint()*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X.adjoint() ;
    } else if ( o == 4){
      T = X.transpose()*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X ;
      T = X.transpose()*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X ;
      T = X.transpose()*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X ;
      T = X.transpose()*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X ;
      }
    }

  T.resize( 0, 0) ;

  return ;

} ;

void transform( int o, Eigen::Ref<Eigen::MatrixXcd> X, Eigen::Ref<Eigen::MatrixXcd> M) {
/*
  Given a tranformation matrix, convert the passed array to or from the 
  basis transformation passed.  For now we will assume the passed matrices
  are square.  

  options :
    o = 0 - Perform X^{t}*M*X
    o = 1 - Perform X*M*X^{t}
    o = 2 - Perform X*M*X
    o = 3 - Perform X^{t}*M*X^{t}
    o = 4 - Perform X^{T}*M*X
*/

  Eigen::MatrixXcd::Index xcol = X.cols(), mcol = M.cols() ;
  Eigen::MatrixXcd T ( xcol, xcol) ;
  
  if ( xcol == mcol){
    if ( o == 0){
      T = X.adjoint()*M ;
      M = T*X ;
    } else if ( o == 1){
      T = X*M ;
      M = T*X.adjoint() ;
    } else if ( o == 2){
      T = X*M ;
      M = T*X ;
    } else if ( o == 3){
      T = X.adjoint()*M ;
      M = T*X.adjoint() ;
    } else if ( o == 4){
      T = X.transpose()*M ;
      M = T*X ;
      }
  } else {
    if ( o == 0){
      T = X.adjoint()*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X ;
      T = X.adjoint()*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X ;
      T = X.adjoint()*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X ;
      T = X.adjoint()*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X ;
    } else if ( o == 1){
      T = X*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X.adjoint() ;
      T = X*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X.adjoint() ;
      T = X*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X.adjoint() ;
      T = X*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X.adjoint() ;
    } else if ( o == 2){
      T = X*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X ;
      T = X*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X ;
      T = X*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X ;
      T = X*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X ;
    } else if ( o == 3){
      T = X.adjoint()*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X.adjoint() ;
      T = X.adjoint()*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X.adjoint() ;
      T = X.adjoint()*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X.adjoint() ;
      T = X.adjoint()*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X.adjoint() ;
    } else if ( o == 4){
      T = X.transpose()*M.block( 0, 0, xcol, xcol) ;
      M.block( 0, 0, xcol, xcol) = T*X ;
      T = X.transpose()*M.block( 0, xcol, xcol, xcol) ;
      M.block( 0, xcol, xcol, xcol) = T*X ;
      T = X.transpose()*M.block( xcol, 0, xcol, xcol) ;
      M.block( xcol, 0, xcol, xcol) = T*X ;
      T = X.transpose()*M.block( xcol, xcol, xcol, xcol) ;
      M.block( xcol, xcol, xcol, xcol) = T*X ;
      }
    }

  T.resize( 0, 0) ;

  return ;

} ;

void oao( int nbasis, std::vector<tei>& iarr, std::vector<tei>& ioarr, Eigen::MatrixXd shf) {
  
/* 
    Given an overlap and two electron integrals.  Convert them to an
    orthogonal basis.
*/
  
    int l_lim ;
    int bra ;
    int ket ;
    int ij ;
    int kl ;
    int oi ;
    int oj ;
    int ok ;
    int ol ;
    int t ;
    tei tmp_tei ;
    double val ;
    double shfprd ;
    double oao_tei ;
    time_dbg oao_time = time_dbg("oao_two_electron_int") ;
  
    /* These four loops generate all unique two electron integrals */
    for ( int i=1; i < nbasis+1; i++) {
      for ( int j=1; j < i+1; j++) {
        for ( int k=1; k < i+1; k++) {
          if ( k == i ){
            l_lim = j + 1 ;
          } else {
            l_lim = k + 1 ;
          }
  
          for ( int l=1; l < l_lim; l++) {
  
            oao_tei = d0 ;
  
            /* The transformation is a four index quantity. */
            for ( int mu=1; mu < nbasis+1; mu++ ) {
              for ( int nu=1; nu < nbasis+1; nu++ ) {
                for ( int lm=1; lm < nbasis+1; lm++ ) {
                  for ( int sg=1; sg < nbasis+1; sg++ ) {
  
                    /* Find the permutation to the unique integral */
                    bra = std::max( mu, nu) ;
                    ket = std::max( lm, sg) ;
                    if ( bra > ket ) {
                      oi = std::max( mu, nu) ;
                      oj = std::min( mu, nu) ;
                      ok = std::max( lm, sg) ;
                      ol = std::min( lm, sg) ;
                    } else if ( ket > bra ) {
                      oi = std::max( lm, sg) ;
                      oj = std::min( lm, sg) ;
                      ok = std::max( mu, nu) ;
                      ol = std::min( mu, nu) ;
                    } else if ( ket == bra ) {
                      oi = bra ;
                      oj = std::max( std::min( mu, nu), std::min( lm, sg)) ;
                      ok = bra ;
                      ol = std::min( std::min( mu, nu), std::min( lm, sg)) ;
                    }
  
                      /* Convert the index into lower triangle of lower triangle */
                      ij = oi*( oi - 1)/2 + oj ;
                      kl = ok*( ok - 1)/2 + ol ;
                      t = ij*( ij-1)/2 + kl ;
                      val = iarr[t-1].r_v() ;
                      shfprd = shf.adjoint()(i-1,mu-1)*shf(nu-1,j-1)*shf.adjoint()(k-1,lm-1)*shf(sg-1,l-1) ;
                      oao_tei += shfprd*val ;
                    } /* End sg */
                  } /* End lm */
                } /* End nu */
              } /* End mu */
  
              /* Save our new oao integrals */
              tmp_tei.set( i-1, j-1, k-1, l-1, oao_tei) ;
              ioarr.push_back(tmp_tei) ;
            } /* End l */
          } /* End k */
        } /* End j */
      } /* End i */
  
    oao_time.end() ;

    return ;
  
    } ;

//  void eulrgrd ( int n_psi, int n_thet, int n_phi, std::vector<double>& w_psi, std::vector<double>& w_thet, 
/*       std::vector<double>& w_phi, std::vector<double>& x_psi, std::vector<double>& x_thet, 
       std::vector<double>& x_phi, int SG){
  
   Set up a grid to integrate over the Euler angles.  The grid and weights are set up to satisfy this equation
   
    If the integral is over O+(3) the euation is
    1 = 1/(8 pi^2) \int_{0}^{2 pi} dpsi \int_{0}^{pi} sin(theta)dtheta int_{0}^{2 pi} dphi
   
    If the integral is over SU(2) the euation is
    1 = 1/(8 pi^2) \int_{0}^{2 pi} dpsi \int_{0}^{pi} sin(theta)dtheta int_{0}^{4 pi} dphi
   
    Where the rotation operator is given by
   
    R( psi, theta, phi) = Exp[-psi S_{z}]*Exp[-theta S_{y}]*Exp[-phi S_{z}]
    R( alpha, beta, gamma) = Exp[-alpha S_{z}]*Exp[-beta S_{y}]*Exp[-gamma S_{z}]
   
  
    nst double tpi=2.0*pi ;
    double phi_lim ;
    double fjunk ;
  
    for ( int i=0; i< n_psi; i++ ){
      w_psi.push_back(d0) ;
      x_psi.push_back(d0) ;
    }
  
    gauleg( d0, tpi, n_psi, x_psi, w_psi) ;
  
    for ( int i=0; i < n_psi; i++ ){
      w_psi[i] = w_psi[i]/tpi ;
    }
   
    for ( int i=0; i< n_thet; i++ ){
      w_thet.push_back(d0) ;
      x_thet.push_back(d0) ;
    }
   
    gauleg( d0, pi, n_thet, x_thet, w_thet) ;
  
    for ( int i=0; i< n_thet; i++ ){
      w_thet[i] = std::sin(x_thet[i])*w_thet[i]/2.0 ;
    }
  
    for ( int i=0; i< n_phi; i++ ){
      w_phi.push_back(d0) ;
      x_phi.push_back(d0) ;
    }
  
    if ( SG == 2 ){
      // Do SU(2) integration 
      phi_lim = tpi*2.0 ;
    } else if ( SG == 3 ) {
      // Do O+(3) integration 
      phi_lim = tpi ;
    }
  
    gauleg( d0, phi_lim, n_psi, x_phi, w_phi) ;
  
    for ( int i=0; i< n_phi; i++ ){
      w_phi[i] = w_phi[i]/phi_lim ;
    }
  
    return ;
  
  } ; */
       
//  void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mo, int nb ){
//  /* 
//   * Apply the complex conjugation operator to the molecular orbitals. 
//   * Return the determinant to be used in the calling routine.
//   * 
//   * It is assumed mo has been dimensioned correctly already.
//   * */
//  
//    Eigen::MatrixXcd tmp ;
//  
//    tmp.resize( 2*nb, 2*nb) ;
//  
//    a.get_mos( tmp) ;
//  
//    mo = tmp.conjugate() ;
//  
//    tmp.resize( 0, 0) ;
//  
//    return  ;
//  
//  } ;
//  
//  void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mo, int nb ){
//  /* 
//   * Apply the complex conjugation operator to the molecular orbitals. 
//   * Return the determinant to be used in the calling routine.
//   * 
//   * It is assumed mo has been dimensioned correctly already.
//   * */
//  
//    double Nnty ;
//  
//    Nnty = pi/2.0e0 ;
//  
//    R_s ( nb, a, mo, d0, Nnty, d0) ;
//  
//    return  ;
//  
//  } ;
//  
//  void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mo, int nb ){
//  /* 
//   * Apply the complex conjugation operator to the molecular orbitals. 
//   * Return the determinant to be used in the calling routine.
//   * 
//   * It is assumed mo has been dimensioned correctly already.
//   * */
//  
//    Eigen::MatrixXcd tmp ;
//    double Nnty ;
//  
//    tmp.resize( 2*nb, 2*nb) ;
//    Nnty = pi/2.0e0 ;
//  
//    a.get_mos( mo) ;
//    tmp = mo.conjugate() ;
//    R_s ( nb, tmp, mo, d0, Nnty, d0) ;
//    tmp.resize( 0, 0) ;
//  
//    return  ;
//  
//  } ;

double pfaffian_A( Eigen::Ref<Eigen::MatrixXd> A) {

/*
  Compute the pfaffian using the Aitken block diagonalization ( LDL^{T}) algorithm.
  CAUTION :: This does not use pivoting to enure a numerically stable
  calculation.  This means that the matrix that is passed in must be well
  conditioned until I can write the algorithm to include pivoting.  Use
  the Householder implementation until then.
*/

  int r, i ;
  double pf = d1 ;
  Eigen::MatrixXd::Index n = A.rows() ;
  Eigen::MatrixXd L ( n, n) ;
  Eigen::MatrixXd P ( n, n) ;
  Eigen::MatrixXd S ( 2, 2) ;

  L.setIdentity() ;

  for ( i = 0; i < n ; i+=2) {
    r = n - i - 2 ;
    S = A.block( i, i, 2, 2) ;
    pf *= S( i, i+1) ;
    if ( i == n - 2){ break ;}
    P.setZero() ;
    P.block( 0, 0, r, r) = Eigen::MatrixXd::Identity ( r, r) ;
    P.block( r, r, i + 2, i + 2) = Eigen::MatrixXd::Identity ( i + 2, i + 2) ;
    P.block( i+2, i, r, 2) = A.block( i+2, i, r, 2)*S.inverse() ;
    L = L*P ;
    A.block( i+2, i+2, r, r) = A.block( i+2, i+2, r, r) + A.block( i+2, i, r, 2)*S.inverse()*A.block( i+2, i, r, 2).transpose() ;
    }

  pf = L.determinant()*pf ;

  return pf ;

} ;

cd pfaffian_A( Eigen::Ref<Eigen::MatrixXcd> A) {

/*
  Compute the pfaffian using the Aitken block diagonalization algorithm.
  CAUTION :: This does not use pivoting to enure a numerically stable
  calculation.  This means that the matrix that is passed in must be well
  conditioned until I can write the algorithm to include pivoting.  Use
  the Householder implementation until then.
*/

  int r, i ;
  cd pf = z1 ;
  Eigen::MatrixXcd::Index n = A.rows() ;
  Eigen::MatrixXcd L ( n, n) ;
  Eigen::MatrixXcd P ( n, n) ;
  Eigen::MatrixXcd S ( 2, 2) ;

  L.setIdentity() ;

  for ( i = 0; i < n ; i+=2) {
    r = n - i - 2 ;
    S = A.block( i, i, 2, 2) ;
    pf *= S( 0, 1) ;
    if ( i == (n - 2)){ break ;}
    P.setZero() ;
    P.block( 0, 0, r, r) = Eigen::MatrixXcd::Identity ( r, r) ;
    P.block( r, r, i + 2, i + 2) = Eigen::MatrixXcd::Identity ( i + 2, i + 2) ;
    P.block( i+2, i, r, 2) = A.block( i+2, i, r, 2)*S.inverse() ;
    L = L*P ;
    A.block( i+2, i+2, r, r) = A.block( i+2, i+2, r, r) + A.block( i+2, i, r, 2)*S.inverse()*A.block( i+2, i, r, 2).transpose() ;
    }

  pf = L.determinant()*pf ;

  return pf ;

} ;

double pfaffian_H( Eigen::Ref<Eigen::MatrixXd> A) {

  int i ;
  double pf = d1, alp, r ;
  Eigen::MatrixXd::Index n = A.rows() ;
  Eigen::MatrixXd P ( n, n) ;
  Eigen::VectorXd x ( n) ;

  x.setZero() ;

  if ( n % 2 == 1){
    return d0 ;
    }

  for ( i = 1; i < n-1; i+=2){
    x.tail( n - i) = A.col( i-1).tail( n - i) ;
    alp = -(x.tail( n - i)).norm() ;
    r = std::sqrt( ( std::pow( alp, d2) - A( i, i-1)*alp)/d2) ;
    x( i) -= alp ;
    x /= d2*r ;
    P.setIdentity() ;
    P -= d2*x*x.adjoint() ;
    A = P*A*P.transpose() ;
    x( i) = d0 ;
    x( i+1) = d0 ;
    pf *= A( i-1, i) ;
    }

  pf *= A( n-2, n-1)*std::pow( -d1, (n-2)/2 ) ;

  x.resize( 0) ;
  P.resize( 0, 0) ;

  return pf ;

} ;

cd pfaffian_H( Eigen::Ref<Eigen::MatrixXcd> A) {
/*
  "Numeric and symbolic evaluation of the pfaffian of general skew-symmetric matrices "
  C. Gonzalez-Ballestero; l.M. Robledo; G. F. Bertsch ;
  Comp. Phys. Comm. 188, 10 2011
  doi:10.1016/j.cpc.2011.04.025

  Numerical Analysis, Burden and Faires, 8th edition
*/
  int i ;
  cd pf = z1, alp, r ;
  Eigen::MatrixXcd::Index n = A.rows() ;
  Eigen::MatrixXcd P ( n, n) ;
  Eigen::VectorXcd x ( n) ;

  x.setZero() ;

  if ( n % 2 == 1){
    return z0 ;
    }

  for ( i = 1; i < n-1; i+=2){
    x.tail( n - i) = A.col( i-1).tail( n - i) ;
    alp = -std::exp( cd( 0.0e0, std::arg( A(i, i-1))))*(x.tail( n - i)).norm() ;
    r = std::sqrt( ( std::pow( alp, d2) - A( i, i-1)*alp)/z2) ;
    x( i) -= alp ;
    x /= z2*r ;
    P.setIdentity() ;
    P -= z2*x*x.adjoint() ;
    A = P*A*P.transpose() ;
    x( i) = z0 ;
    x( i+1) = z0 ;
    pf *= A( i-1, i) ;
    }

  pf *= A( n-2, n-1)*std::pow( -z1, (n-2)/2 ) ;

  x.resize( 0) ;
  P.resize( 0, 0) ;

  return pf ;

} ;

void rand_unitary( Eigen::Ref<Eigen::MatrixXcd> u) {
/*
  Generate a random unitary matrix with Haar distribution
*/
  int ur, uc ;
  double rl, im ;
  unsigned long int seed = 7652413 ;
  time_t clock ;
  cd val ;
  const gsl_rng_type* T ;
  gsl_rng* r ;
  seed *= static_cast<unsigned long int>(time(&clock) % 100) ;
  Eigen::MatrixXcd::Index cols = u.cols(), rows = u.rows() ;
  Eigen::VectorXcd v( rows) ;
  Eigen::MatrixXcd Q( rows, cols), R( rows, cols), L( rows, cols) ;
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr ;

  gsl_rng_env_setup() ;

  T = gsl_rng_ranlxs0 ;
  r = gsl_rng_alloc( T) ;
  gsl_rng_set( r, seed) ;

  for( ur = 0; ur < rows; ur++){
    for( uc = 0; uc < cols; uc++){
      rl = d2*gsl_rng_uniform_pos(r) - d1 ;
      im = d2*gsl_rng_uniform_pos(r) - d1 ;
      val = std::complex<double> ( rl, im) ;
      u( ur, uc) = val ;
      }
    }

  qr.compute( u) ;

  Q = qr.householderQ() ;
  R = qr.matrixQR() ;
  v = R.diagonal() ;
  for ( ur = 0; ur < rows; ur++){
    v(ur) = v(ur)/abs(v(ur)) ;
    }
  L = v.asDiagonal() ;

  u = Q*L ;

  gsl_rng_free( r) ;
  Q.resize( 0, 0) ;
  R.resize( 0, 0) ;
  L.resize( 0, 0) ;
  v.resize( 0) ;
  

  return ;

  } ;

void testing_magnetic_structure( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> mx, Eigen::Ref<Eigen::MatrixXcd> my, Eigen::Ref<Eigen::MatrixXcd> mz){

/*
  Given the magnetization densities determine whether the wavefunction is collinear or coplanar
*/

  int i, mchk, nzero, n = mx.cols() ;
  double thresh = 1.0e-6 ;
  bool preal = false, mxreal = false, myreal = false, mzreal = false ;
  bool mxzero = false, myzero = false, mzzero = false ;
  bool colin = false, ncopla = false ;
  cd test ;
  Eigen::Matrix3d Tr ;
  Eigen::Matrix3cd T ;
  Eigen::Matrix3cd evec ;
  Eigen::Vector3cd v ;
  Eigen::MatrixXd scr3( n, n) ;
  Eigen::MatrixXcd scr( n, n), scr1( n, n), scr2( n, n) ;
  Eigen::EigenSolver<Eigen::Matrix3d> Tr_diag ;
  Eigen::ComplexEigenSolver<Eigen::Matrix3cd> T_diag ;

  scr = mx*mx ;
  T( 0, 0) = scr.trace() ;
  scr = my*mx ;
  T( 1, 0) = scr.trace() ;
  scr = mz*mx ;
  T( 2, 0) = scr.trace() ;
  scr = mx*my ;
  T( 0, 1) = scr.trace() ;
  scr = mx*mz ;
  T( 0, 2) = scr.trace() ;
  scr = my*my ;
  T( 1, 1) = scr.trace() ;
  scr = mz*my ;
  T( 2, 1) = scr.trace() ;
  scr = my*mz ;
  T( 1, 2) = scr.trace() ;
  scr = mz*mz ;
  T( 2, 2) = scr.trace() ;

  T_diag.compute( T) ;
  v = T_diag.eigenvalues() ;
  evec = T_diag.eigenvectors() ;

/*
  Rotate the magnetization densities so that 
  ||my|| < ||mx|| < ||mz||
*/
  scr  = evec( 0, 0)*mx + evec( 0, 1)*my + evec( 0, 2)*mz ;
  scr1 = evec( 1, 0)*mx + evec( 1, 1)*my + evec( 1, 2)*mz ;
  scr2 = evec( 2, 0)*mx + evec( 2, 1)*my + evec( 2, 2)*mz ;

  mx = scr1 ;
  my = -scr ;
  mz = scr2 ;

/*
  Determine whether they are zero
*/
  scr3 = mx.cwiseAbs() ;
  mxzero = ( scr3.array() < 1.0e-8).all() ;
  scr3 = my.cwiseAbs() ;
  myzero = ( scr3.array() < 1.0e-8).all() ;
  scr3 = mz.cwiseAbs() ;
  mzzero = ( scr3.array() < 1.0e-8).all() ;
/*
  Determine whether they are real of complex matrices
*/
  scr = p - p.conjugate() ;
  scr3 = scr.imag().cwiseAbs() ;
  preal = ( scr3.array() < 1.0e-8).all() ;

  if ( ! mxzero ){
    scr = mx - mx.conjugate() ;
    scr3 = scr.imag().cwiseAbs() ;
    mxreal = ( scr3.array() < 1.0e-8).all() ;
    }

  if ( ! myzero ){
    scr = my - my.conjugate() ;
    scr3 = scr.imag().cwiseAbs() ;
    myreal = ( scr3.array() < 1.0e-8).all() ;
    }

  if ( ! mzzero ){
    scr = mz - mz.conjugate() ;
    scr3 = scr.imag().cwiseAbs() ;
    mzreal = ( scr3.array() < 1.0e-8).all() ;
    }

/*
  Determine the number of zeros
*/

  nzero= 0 ;
  for( i = 0; i < 3; i++){
    test = std::abs(v(i)) ;
    if ( std::real( test) < thresh){
      nzero ++ ;
      mchk = static_cast<int>(std::floor(std::log10(std::real( test)))) + 10 ;
      if ( mchk == 1){
        std::cout << " Caution: Eigenvalue is below the thresh hold by no more than a factor of 10 " << std::endl ;
        }
      }
    std::cout << "       eig           abs(eig) : " << std::endl ;
    std::cout << v(i) << "    " << test << std::endl ;
    }

  if ( nzero == 3){
     ;
  } else if ( nzero == 2){
    colin = true ;
  } else if ( nzero < 2){
/*
  Let's check for coplanarity since we have three non-zero eigenvalues
*/
    scr.setZero() ;
    scr.real() = mz.real()*mz.real() ;
    Tr( 0, 0) = scr.real().trace() ;
    scr.real() = my.real()*mz.real() ;
    Tr( 1, 0) = scr.real().trace() ;
    scr.real() = mx.real()*mz.real() ;
    Tr( 2, 0) = scr.real().trace() ;
    scr.real() = mz.real()*my.real() ;
    Tr( 0, 1) = scr.real().trace() ;
    scr.real() = mz.real()*mx.real() ;
    Tr( 0, 2) = scr.real().trace() ;
    scr.real() = my.real()*my.real() ;
    Tr( 1, 1) = scr.real().trace() ;
    scr.real() = mx.real()*my.real() ;
    Tr( 2, 1) = scr.real().trace() ;
    scr.real() = my.real()*mx.real() ;
    Tr( 1, 2) = scr.real().trace() ;
    scr.real() = mx.real()*mx.real() ;
    Tr( 2, 2) = scr.real().trace() ;

    Tr_diag.compute( Tr) ;
    v = Tr_diag.eigenvalues() ;

    /*
      Generate new magnetization densities
    */

    nzero = 0 ;
    for( i = 0; i < 3; i++){
      test = std::abs(v(i)) ;
      if ( std::real( test) < thresh){
        nzero++ ;
        mchk = static_cast<int>(std::floor(std::log10(std::real( test)))) + 10 ;
        if ( mchk == 1) {
          std::cout << " Caution: Eigenvalue is below the thresh hold by no more than a factor of 10 " << std::endl ;
          }
        }
      std::cout << "       eig           abs(eig) : " << std::endl ;
      std::cout << v(i) << "    " << test << std::endl ;
      }
    }
    if ( nzero == 0){
      ncopla = true ;
      }
  /*
    Determine the character of the wavefunction now that we have all the information we need.
  */
  std::cout << " Characterizing the Wavefunction " << std::endl ;
  std::cout << " -                             - " << std::endl ;
  if ( mxzero && myzero && mzzero && preal ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   TICS        real RHF       S^{2},S_{z},K,O" << std::endl ;
  } else if ( mxzero && myzero && mzzero && ! preal ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   CCW       complex RHF       S^{2},S_{z}" << std::endl ;
  } else if ( mxzero && myzero && ! mzreal && preal ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   ASCW      paired UHF           S_{z},O" << std::endl ;
  } else if ( mxzero && myzero && mzreal && preal ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   ASDW        real UHF           S_{z},K" << std::endl ;
  } else if ( mxzero && myzero ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   ASW       complex UHF           S_{z}" << std::endl ;
  } else if ( ! mxreal && ! myreal && ! mzreal && preal ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   TSCW      paired GHF             O" << std::endl ;
  } else if ( mxreal && myreal && mzreal && preal ){
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   TSDW       real GHF              K" << std::endl ;
  } else {
    std::cout << " Fukutome   Stuber-Paldus       Symmetries" << std::endl ;
    std::cout << "   TSW       complex GHF             " << std::endl ;
    if ( colin){
      std::cout << "   Wavefunction is collinear " << std::endl ;
    } else if ( ncopla){
      std::cout << "   Wavefunction is non-coplanar " << std::endl ;
    } else {
      std::cout << "   Wavefunction is coplanar " << std::endl ;
      }
    }

  return ;

  } ;

void clean_mat( Eigen::Ref<Eigen::MatrixXcd> m){
/*
  Given a matrix "clean it" by zeroing small numbers
  Currently only works for complex double.
*/
  int i ;
  double tr, ti ;
  cd t ;
  for ( i = 0; i < m.size(); i++){
    t = *(m.data() + i) ;
    tr = t.real() ;
    ti = t.imag() ;
    if ( tr < 1.0e-13){
      tr = d0 ;
      }
    if ( ti < 1.0e-13){
      ti = d0 ;
      }
    *(m.data() + i) = cd ( tr, ti) ;
    }

  return ;

} ;

uint32_t xorshift32( uint32_t *state){
  uint32_t x = *state ;

  x ^= x << 13 ;
  x ^= x >> 17 ;
  x ^= x << 5 ;
  *state = x ;

  return x ;

  } ;

double rand01( uint32_t *state) {

/*
  A pseudo-random number generator that can be quickly implemented in 
  multiple environments for debugging "randomly" initialized matrices.
*/

  uint32_t v ;
  v = xorshift32( state) ;
  return static_cast<double>( v)/static_cast<double>(UINT32_MAX) - 0.5e0 ;

  } ;

void wavefunction_characterization( int wt, const int& nbas, const int& nalp, const int& nbet, const std::string& f) {
/*
  Given the name of binary file open it for testing the wavefunction.
  If no name is given then open the default.

  Note : this is only used for GHF type wavefunctions and so the wavefunction
  that is read in must be GHF.

  wt - what type of wavefunction are we loading?
*/
  int nele = nalp + nbet ;
  Eigen::MatrixXcd p( 2*nbas, 2*nbas), c( nbas, nbas), mx( nbas, nbas), my( nbas, nbas), mz( nbas, nbas) ;
  p.setZero() ;

  if ( wt % 2 == 1) {
    wfn< double, Eigen::Dynamic, Eigen::Dynamic> w ;
    if ( (wt + 1)/2 == 1 ) {
      w.eig.resize( nbas) ;
      w.moc.resize( nbas, nbas) ;
      load_wfn( w) ;
      p.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nalp)*w.moc.block( 0, 0, nbas, nalp).adjoint() ;
      p.block( nbas, nbas, nbas, nbas) = w.moc.block( 0, 0, nbas, nalp)*w.moc.block( 0, 0, nbas, nalp).adjoint() ;
    } else if ( (wt + 1)/2 == 2 ) {
      w.eig.resize( 2*nbas) ;
      w.moc.resize( nbas, 2*nbas) ;
      load_wfn( w) ;
      p.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nalp)*w.moc.block( 0, 0, nbas, nalp).adjoint() ;
      p.block( nbas, nbas, nbas, nbas) = w.moc.block( 0, nbas, nbas, nbet)*w.moc.block( 0, nbas, nbas, nbet).adjoint() ;
    } else {
      w.eig.resize( 2*nbas) ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      load_wfn( w) ;
      p = w.moc.block( 0, 0, 2*nbas, nele)*w.moc.block( 0, 0, 2*nbas, nele).adjoint() ;
      }
  } else {
    wfn< cd, Eigen::Dynamic, Eigen::Dynamic> w ;
    if ( (wt + 1)/2 == 1 ) {
      w.eig.resize( nbas) ;
      w.moc.resize( nbas, nbas) ;
      load_wfn( w) ;
      p.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nalp)*w.moc.block( 0, 0, nbas, nalp).adjoint() ;
      p.block( nbas, nbas, nbas, nbas) = w.moc.block( 0, 0, nbas, nalp)*w.moc.block( 0, 0, nbas, nalp).adjoint() ;
    } else if ( (wt + 1)/2 == 2 ) {
      w.eig.resize( 2*nbas) ;
      w.moc.resize( nbas, 2*nbas) ;
      load_wfn( w) ;
      p.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nalp)*w.moc.block( 0, 0, nbas, nalp).adjoint() ;
      p.block( nbas, nbas, nbas, nbas) = w.moc.block( 0, nbas, nbas, nbet)*w.moc.block( 0, nbas, nbas, nbet).adjoint() ;
    } else {
      w.eig.resize( 2*nbas) ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      load_wfn( w) ;
      p = w.moc.block( 0, 0, 2*nbas, nele)*w.moc.block( 0, 0, 2*nbas, nele).adjoint() ;
      }
    }
  
  charge_magnetic_decomposition( p, c, mx, my, mz) ;

  testing_magnetic_structure( p, mx, my, mz) ;

  return ;

  } ;

