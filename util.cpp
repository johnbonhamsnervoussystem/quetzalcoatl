#include "constants.h"
#include<Eigen/Dense>
#include<Eigen/LU>
#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include "hfwfn.h"
#include "solver.h"
#include "util.h"
#include "integr.h"
#include "tei.h"
#include "wigner.h"

/* Utilities that don't belong elsewhere. */

float fact ( int n) {
  /* Return the factorial of n */
  float f=1.0 ;

  if ( n < 1 ) { return 1.0 ;}

  while ( n > 1 ) {
    f = f*float(n) ;
    n+= -1 ;
  }

  return f ;

} ;

void oao( int nbasis, hfwfn& a, Eigen::MatrixXf s) {

  /* Given an overlap and a Slater determinant, check if the determinant is
   * already in the orthogonal atomic orbital basis.  If not, put it into the 
   * oao basis. */

  int iwfnt ;
  float detr ;
  float detru ;
  float d1 = 1e0 ;
  float thresh = 1e-7 ;
  std::complex<float> detc ;
  std::complex<float> detcu ;
  Eigen::MatrixXf mor ; 
  Eigen::MatrixXf tmpr ; 
  Eigen::MatrixXf shf_r ; 
  Eigen::MatrixXcf moc ; 
  Eigen::MatrixXcf tmpc ; 
  Eigen::MatrixXcf shf_c ; 
  
  iwfnt = a.get_wti() ;
  if ( iwfnt % 2 == 1 ){

  /* Real wavefuntion*/

    if ( iwfnt == 1 ){

    /* Restricted wavefuntion*/

      mor.resize( nbasis, nbasis) ;
      tmpr.resize( nbasis, nbasis) ;
      a.get_mos( mor) ;
      tmpr = mor.adjoint()*mor ;
      detr = tmpr.determinant() ;

      if ( std::abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr = shf_r.adjoint()*s*mor ;

        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

    } else if ( iwfnt == 3 ){

    /* unrestricted wavefuntion*/

      mor.resize( nbasis, 2*nbasis) ;
      tmpr.resize( nbasis, 2*nbasis) ;
      a.get_mos( mor) ;
      tmpr.block( 0, 0, nbasis, nbasis) = mor.block( 0, 0, nbasis, nbasis).adjoint()*mor.block( 0, 0, nbasis, nbasis) ;
      tmpr.block( 0, nbasis, nbasis, nbasis) = mor.block( 0, nbasis, nbasis, nbasis).adjoint()*mor.block( 0, nbasis, nbasis, nbasis) ;
      detr = tmpr.block( 0, 0, nbasis, nbasis).determinant() ;
      detru = tmpr.block( 0, nbasis, nbasis, nbasis).determinant() ;

      if ( std::abs(detr - d1) > thresh |  std::abs(detru - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr.block( 0, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, 0, nbasis, nbasis) ;
        tmpr.block( 0, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

    } else if ( iwfnt == 5 ){

    /* Generalized wavefuntion*/

      mor.resize( 2*nbasis, 2*nbasis) ;
      tmpr.resize( 2*nbasis, 2*nbasis) ;
      a.get_mos( mor) ;
      tmpr = mor.adjoint()*mor ;
      detr = tmpr.determinant() ;

      if ( std::abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr.block( 0, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, 0, nbasis, nbasis) ;
        tmpr.block( 0, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, nbasis, nbasis, nbasis) ;
        tmpr.block( nbasis, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( nbasis, 0, nbasis, nbasis) ;
        tmpr.block( nbasis, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( nbasis, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;


      }

    } else {

    /* complex wavefuntion*/

    if ( iwfnt == 2 ){

    /* Restricted wavefuntion*/

      moc.resize( nbasis, nbasis) ;
      tmpc.resize( nbasis, nbasis) ;
      a.get_mos( moc) ;
      tmpc = moc.adjoint()*moc ;
      detc = tmpc.determinant() ;
      detr = (std::conj( detc)*detc).real() ;

      if ( std::abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        tmpc = shf_c.adjoint()*s*mor ;
        
        a.set_mos( tmpc) ;

      } 

      shf_c.resize( 0, 0) ; 
      tmpc.resize( 0, 0) ;
      moc.resize( 0, 0) ;

    } else if ( iwfnt == 4 ){

    /* unrestricted wavefuntion*/

      moc.resize( nbasis, 2*nbasis) ;
      tmpc.resize( nbasis, 2*nbasis) ;
      a.get_mos( moc) ;
      tmpc.block( 0, 0, nbasis, nbasis) = moc.block( 0, 0, nbasis, nbasis).adjoint()*moc.block( 0, 0, nbasis, nbasis) ;
      tmpc.block( 0, nbasis, nbasis, nbasis) = moc.block( 0, nbasis, nbasis, nbasis).adjoint()*moc.block( 0, nbasis, nbasis, nbasis) ;
      detc = tmpc.block( 0, 0, nbasis, nbasis).determinant() ;
      detcu = tmpc.block( 0, nbasis, nbasis, nbasis).determinant() ;
      detr = (std::conj( detc)*detc).real() ;
      detru = (std::conj( detcu)*detcu).real() ;

      if ( std::abs(detr - d1) > thresh |  std::abs(detru - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        tmpc.block( 0, 0, nbasis, nbasis) = shf_c.adjoint()*s*moc.block( 0, 0, nbasis, nbasis) ;
        tmpc.block( 0, nbasis, nbasis, nbasis) = shf_c.adjoint()*s*moc.block( 0, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpc) ;

      } 

      shf_c.resize( 0, 0) ; 
      tmpc.resize( 0, 0) ;
      moc.resize( 0, 0) ;

    } else if ( iwfnt == 6 ){

    /* Generalized wavefuntion*/

      moc.resize( 2*nbasis, 2*nbasis) ;
      tmpc.resize( 2*nbasis, 2*nbasis) ;
      a.get_mos( moc) ;
      tmpc = moc.adjoint()*moc ;
      detc = tmpc.determinant() ;
      detr = (std::conj( detc)*detc).real() ;

      if ( std::abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        tmpc.block( 0, 0, nbasis, nbasis) = shf_c.adjoint()*s*moc.block( 0, 0, nbasis, nbasis) ;
        tmpc.block( 0, nbasis, nbasis, nbasis) = shf_c.adjoint()*s*moc.block( 0, nbasis, nbasis, nbasis) ;
        tmpc.block( nbasis, 0, nbasis, nbasis) = shf_c.adjoint()*s*moc.block( nbasis, 0, nbasis, nbasis) ;
        tmpc.block( nbasis, nbasis, nbasis, nbasis) = shf_c.adjoint()*s*moc.block( nbasis, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpc) ;

      } 

      shf_c.resize( 0, 0) ; 
      tmpc.resize( 0, 0) ;
      moc.resize( 0, 0) ;

      }  

    }

  }  ;

void oao( int nbasis, Eigen::Ref<Eigen::MatrixXf> ouv, Eigen::MatrixXf s) {

  /* Given an overlap and matrix elements, put the matrix elements into
  *  the orthogonal ao basis.  Since our basis will alwyas be real, the
  *  routine will only accept the nbasis*nbasis real matrix.  Any conversion
  *  to complex or ghf must occur in the calling routines.*/

  Eigen::MatrixXf tmpr ;
  Eigen::MatrixXf shf_r ;
  Eigen::MatrixXcf shf_c ;

  shf_r.resize( nbasis, nbasis) ;
  shf_c.resize( nbasis, nbasis) ;
  tmpr.resize( nbasis, nbasis) ;

  canort( s, shf_c, nbasis) ;
  shf_r = shf_c.real() ;

  tmpr = ouv*shf_r ;
  ouv = shf_r.adjoint()*tmpr ;

  tmpr.resize( 0, 0) ;
  shf_c.resize( 0, 0) ;
  shf_r.resize( 0, 0) ;

  return ;

  } ;

void oao( int nbasis, std::vector<tei>& iarr, std::vector<tei>& ioarr, Eigen::MatrixXf s) {

  /* 
   * Given an overlap and two electron integrals.  Convert them to an
   * orthogonal basis.
   */

  Eigen::MatrixXf shf_r ;
  Eigen::MatrixXcf shf_c ;
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
  float val ;
  float shfprd ;
  float oao_tei ;

  shf_r.resize( nbasis, nbasis) ;
  shf_c.resize( nbasis, nbasis) ;

  canort( s, shf_c, nbasis) ;
  shf_r = -1*shf_c.real() ;

  shf_c.resize( 0, 0) ;
  /* I imagine there is a better way to do this but you don't know
   * what you don't know. */

  /* These four loops generate all uniqure two electron integrals */
  for ( int i=1; i < nbasis+1; i++) {
    for ( int j=1; j < i+1; j++) {
      for ( int k=1; k < i+1; k++) {
        if ( k == i ){
          l_lim = j + 1 ;
        } else {
          l_lim = k + 1 ;
        }

        for ( int l=1; l < l_lim; l++) {

          oao_tei =  0.0 ;

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
                    shfprd = shf_r.adjoint()(i-1,mu-1)*shf_r(nu-1,j-1)*shf_r.adjoint()(k-1,lm-1)*shf_r(sg-1,l-1) ;
                    oao_tei += shfprd*val ;
                  } /* End sg */
                } /* End lm */
              } /* End nu */
            } /* End mu */

            /* Save our new oao integrals */
            tmp_tei.set( i, j, k, l, oao_tei) ;
            ioarr.push_back(tmp_tei) ;
          } /* End l */
        } /* End k */
      } /* End j */
    } /* End i */

  shf_r.resize( 0, 0) ;

  return ;

  } ;

void eulrgrd ( int n_psi, int n_thet, int n_phi, std::vector<float>& w_psi, std::vector<float>& w_thet, 
     std::vector<float>& w_phi, std::vector<float>& x_psi, std::vector<float>& x_thet, 
     std::vector<float>& x_phi, int SG){

/* Set up a grid to integrate over the Euler angles.  The grid and weights are set up to satisfy this equation
 *
 * If the integral is over O+(3) the euation is
 * 1 = 1/(8 pi^2) \int_{0}^{2 pi} dpsi \int_{0}^{pi} sin(theta)dtheta int_{0}^{2 pi} dphi
 *
 * If the integral is over SU(2) the euation is
 * 1 = 1/(8 pi^2) \int_{0}^{2 pi} dpsi \int_{0}^{pi} sin(theta)dtheta int_{0}^{4 pi} dphi
 *
 * Where the rotation operator is given by
 *
 * R( psi, theta, phi) = Exp[-psi S_{z}]*Exp[-theta S_{y}]*Exp[-phi S_{z}]
 * R( alpha, beta, gamma) = Exp[-alpha S_{z}]*Exp[-beta S_{y}]*Exp[-gamma S_{z}]
 *
 * */

  const float tpi=2.0*pi ;
  float phi_lim ;
  float fjunk ;

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
    /* Do SU(2) integration */
    phi_lim = tpi*2.0 ;
  } else if ( SG == 3 ) {
    /* Do O+(3) integration */
    phi_lim = tpi ;
  }

  gauleg( d0, phi_lim, n_psi, x_phi, w_phi) ;

  for ( int i=0; i< n_phi; i++ ){
    w_phi[i] = w_phi[i]/phi_lim ;
  }

  return ;

} ;

void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcf> mo, int nb ){
/* 
 * Apply the complex conjugation operator to the molecular orbitals. 
 * Return the determinant to be used in the calling routine.
 * 
 * It is assumed mo has been dimensioned correctly already.
 * */

  Eigen::MatrixXcf tmp ;

  tmp.resize( 2*nb, 2*nb) ;

  a.get_mos( tmp) ;

  mo = tmp.conjugate() ;

  tmp.resize( 0, 0) ;

  return  ;

} ;

void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcf> mo, int nb ){
/* 
 * Apply the complex conjugation operator to the molecular orbitals. 
 * Return the determinant to be used in the calling routine.
 * 
 * It is assumed mo has been dimensioned correctly already.
 * */

  float Nnty ;

  Nnty = pi/2.0e0 ;

  R_s ( nb, a, mo, d0, Nnty, d0) ;

  return  ;

} ;

void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcf> mo, int nb ){
/* 
 * Apply the complex conjugation operator to the molecular orbitals. 
 * Return the determinant to be used in the calling routine.
 * 
 * It is assumed mo has been dimensioned correctly already.
 * */

  Eigen::MatrixXcf tmp ;
  float Nnty ;

  tmp.resize( 2*nb, 2*nb) ;
  Nnty = pi/2.0e0 ;

  a.get_mos( mo) ;
  tmp = mo.conjugate() ;
  R_s ( nb, tmp, mo, d0, Nnty, d0) ;
  tmp.resize( 0, 0) ;

  return  ;

} ;

