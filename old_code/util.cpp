#include <algorithm>
#include <cmath>
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include "hfwfn.h"
#include "integr.h"
#include <iostream>
#include "solver.h"
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include <vector>
#include "wigner.h"

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
    n = static_cast<double>(k) + 3.0e0/d2 ;
    f = gsl_sf_hyperg_1F1( m, n, -t)/( d2*static_cast<double>(k) + d1) ;
    
    return f ;

  } ;

//template<class matrix>
//void oao( int& nbasis, int& wfntyp, matrix& a, Eigen::Ref<Eigen::MatrixXd> s, Eigen::Ref<Eigen::MatrixXd> x) {
/* 
  Given a Slater determinant put it into the oao basis.
*/
//    typename matrix::Scalar f ;
//    matrix t ; 
//    time_dbg oao_wfn = time_dbg("oao_wfn") ;
//
//    std::cout << f << std::endl ;
//
//    if ( wfntyp % 2 == 1) {
//      /* Real wavefuntion*/
//      if ( wfntyp == 1) {
//        /* Restricted wavefuntion*/
//        t.resize( nbasis, nbasis) ;
//        t = x.adjoint()*s*a ;
//        a = t ;
//  
//      } else if ( wfntyp == 3 ){
//        /* unrestricted wavefuntion*/
//        t.resize( nbasis, 2*nbasis) ;
//        t.block( 0, 0, nbasis, nbasis) = x.adjoint()*s*a.block( 0, 0, nbasis, nbasis) ;
//        t.block( 0, nbasis, nbasis, nbasis) = x.adjoint()*s*a.block( 0, nbasis, nbasis, nbasis) ;
//        a = t ;
//  
//      } else if ( wfntyp == 5 ){
//        /* Generalized wavefuntion*/
//        t.resize( 2*nbasis, 2*nbasis) ;
//        t.block( 0, 0, nbasis, nbasis) = x.adjoint()*s*a.block( 0, 0, nbasis, nbasis) ;
//        t.block( 0, nbasis, nbasis, nbasis) = x.adjoint()*s*a.block( 0, nbasis, nbasis, nbasis) ;
//        t.block( nbasis, 0, nbasis, nbasis) = x.adjoint()*s*a.block( nbasis, 0, nbasis, nbasis) ;
//        t.block( nbasis, nbasis, nbasis, nbasis) = x.adjoint()*s*a.block( nbasis, nbasis, nbasis, nbasis) ;
//    
//        a = t ;
//        }
//  
//        t.resize( 0, 0) ;
//      } else {
//      /* Real wavefuntion*/
//        if ( wfntyp == 2) {
//          /* Restricted wavefuntion*/
//          t.resize( nbasis, nbasis) ;
//          t = x.adjoint()*s*a ;
//          a = t ;
//  
//        } else if ( wfntyp == 4 ){
//          /* unrestricted wavefuntion*/
//          t.resize( nbasis, 2*nbasis) ;
//          t.block( 0, 0, nbasis, nbasis) = x.adjoint()*s*a.block( 0, 0, nbasis, nbasis) ;
//          t.block( 0, nbasis, nbasis, nbasis) = x.adjoint()*s*a.block( 0, nbasis, nbasis, nbasis) ;
//          a = t ;
//  
//        } else if ( wfntyp == 6 ){
//          /* Generalized wavefuntion*/
//          t.resize( 2*nbasis, 2*nbasis) ;
//          t.block( 0, 0, nbasis, nbasis) = x.adjoint()*s*a.block( 0, 0, nbasis, nbasis) ;
//          t.block( 0, nbasis, nbasis, nbasis) = x.adjoint()*s*a.block( 0, nbasis, nbasis, nbasis) ;
//          t.block( nbasis, 0, nbasis, nbasis) = x.adjoint()*s*a.block( nbasis, 0, nbasis, nbasis) ;
//          t.block( nbasis, nbasis, nbasis, nbasis) = x.adjoint()*s*a.block( nbasis, nbasis, nbasis, nbasis) ;
//    
//          a = t ;
//          }
//  
//          t.resize( 0, 0) ;
//  
//        }
//
//    oao_wfn.end() ;
//
//    return ;
//  
//  }  ;

//template void oao( int&, int&, Eigen::Ref<Eigen::MatrixXd>&, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>) ;
//template void oao( int&, int&, Eigen::Ref<Eigen::MatrixXcd>&, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>) ;

/*void oao( Eigen::Ref<Eigen::MatrixXd> tmp, Eigen::Ref<Eigen::MatrixXd> O, Eigen::MatrixXd X) {
*/
/*
  This is a wrapper to perform a basic function and make the code clearner.
  Pass this two matrices and a transformation and do a similary transform.

*/
/*    tmp = X.adjoint()*O*X ;
    O = tmp ;

    return ;

    } ;*/

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
       
  void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mo, int nb ){
  /* 
   * Apply the complex conjugation operator to the molecular orbitals. 
   * Return the determinant to be used in the calling routine.
   * 
   * It is assumed mo has been dimensioned correctly already.
   * */
  
    Eigen::MatrixXcd tmp ;
  
    tmp.resize( 2*nb, 2*nb) ;
  
    a.get_mos( tmp) ;
  
    mo = tmp.conjugate() ;
  
    tmp.resize( 0, 0) ;
  
    return  ;
  
  } ;
  
  void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mo, int nb ){
  /* 
   * Apply the complex conjugation operator to the molecular orbitals. 
   * Return the determinant to be used in the calling routine.
   * 
   * It is assumed mo has been dimensioned correctly already.
   * */
  
    double Nnty ;
  
    Nnty = pi/2.0e0 ;
  
    R_s ( nb, a, mo, d0, Nnty, d0) ;
  
    return  ;
  
  } ;
  
  void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> mo, int nb ){
  /* 
   * Apply the complex conjugation operator to the molecular orbitals. 
   * Return the determinant to be used in the calling routine.
   * 
   * It is assumed mo has been dimensioned correctly already.
   * */
  
    Eigen::MatrixXcd tmp ;
    double Nnty ;
  
    tmp.resize( 2*nb, 2*nb) ;
    Nnty = pi/2.0e0 ;
  
    a.get_mos( mo) ;
    tmp = mo.conjugate() ;
    R_s ( nb, tmp, mo, d0, Nnty, d0) ;
    tmp.resize( 0, 0) ;
  
    return  ;
  
  } ;

double pfaffian( Eigen::Ref<Eigen::MatrixXd> m){

  /* Currently just double.  Templated version working below. */
  Eigen::MatrixXd::Index i = m.rows() ;
  Eigen::MatrixXd::Index j = m.cols() ;
  int k, h ;
  const size_t N = i ;
  double sgn, pf ;
  double x = 0.0 ;

  /* Do some error handling */

  if ( i != j ){
    std::cout << " Non square matrix passed to pfaffian: " << std::endl ;
    std::cout << " 1.0 returned " << std::endl ;
    return 1.0 ;
    }

  if ( i > GSL_SF_FACT_NMAX ){
    std::cout << " GSL factorial is limited to values <= " << GSL_SF_FACT_NMAX << std::endl ;
    std::cout << " Tell Kyle to add gamma functionality for larger systems. " << std::endl ;
    std::cout << " 1.0 returned " << std::endl ;
    return 1.0 ;
    }

  h = i/2 ;

  gsl_permutation* p = gsl_permutation_alloc(N) ;
  gsl_permutation_init(p) ;

  do {
    sgn = pow( -1.0, gsl_permutation_inversions( p)) ;
    pf = 1.0 ;
    for( k=1; k <= h; k++){
      pf *= m( p->data[2*k-2], p->data[2*k-1]) ;
      }
    x += sgn*pf ;
    }
  while( gsl_permutation_next(p) == GSL_SUCCESS) ;

  gsl_permutation_free(p) ;

  sgn = pow( 2.0, h)*gsl_sf_fact( h) ;
  pf = x/sgn ;

  return pf ;

} ;

cd pfaffian( Eigen::Ref<Eigen::MatrixXcd> m){

  Eigen::MatrixXcd::Index i = m.rows() ;
  Eigen::MatrixXcd::Index j = m.cols() ;
  int k, h ;
  const size_t N = i ;
  cd sgn, pf ;
  cd x = z0 ;

  /* Do some error handling */

  if ( i != j ){
    std::cout << " Non square matrix passed to pfaffian: " << std::endl ;
    std::cout << " 1.0 returned " << std::endl ;
    return 1.0 ;
    }

  if ( i > GSL_SF_FACT_NMAX ){
    std::cout << " GSL factorial is limited to values <= " << GSL_SF_FACT_NMAX << std::endl ;
    std::cout << " Tell Kyle to add gamma functionality for larger systems. " << std::endl ;
    std::cout << " 1.0 returned " << std::endl ;
    return 1.0 ;
    }

  h = i/2 ;

  gsl_permutation* p = gsl_permutation_alloc(N) ;
  gsl_permutation_init(p) ;

  do {
    sgn = z1*pow( -1.0, gsl_permutation_inversions( p)) ;
    pf = z1 ;
    for( k=1; k <= h; k++){
      pf *= m( p->data[2*k-2], p->data[2*k-1]) ;
      }
    x += sgn*pf ;
    }
  while( gsl_permutation_next(p) == GSL_SUCCESS) ;

  gsl_permutation_free(p) ;

  sgn = z1*pow( 2.0, h)*gsl_sf_fact( h) ;
  pf = x/sgn ;

  return pf ;

} ;

