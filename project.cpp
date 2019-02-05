#include <cmath>
#include "common.h"
#include <complex>
#include "constants.h"
#include <deque>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "evalm.h"
#include <iostream>
#include <fstream>
#include "integr.h"
#include "guess.h"
#include "numer.h"
#include "project.h"
#include <Eigen/QR>
#include "qtzio.h"
#include "qtzcntrl.h"
#include "solver.h"
#include "tei.h"
#include <time.h>
#include "time_dbg.h"
#include "util.h"
#include <vector>
#include "diis.h"
#include "wigner.h"
#include "wfn.h"
#include <gsl/gsl_rng.h>

/* Controlling routines to do Projection methods. 

 Wish list -
   -HF
     - Spin
     - C.C.
     - Point Group
   -HFB
     - All of the above
     - Number

 */

void prj_drv( common& com, std::vector<tei>& intarr, int opt){
/*
  Driver routine for projected methods.  This routine will unpack the 
  options and set flags for passing into the drivers for various flavors
  of the model.  

  We need parse a lot of options here.  
    -Type of projection
      - Particle Number
      - Spin
      - Complex Conjugation
      - Time Reversal
    -Number of grid points
      - alpha
      - beta
      - gamma
      - part
    -Type of grid
      - ???
    -Type of wavefunction
      - HF
      - HFB

  I mean, we only got HFB for now!!

  Let's put everything in the orthonormal ao basis here since it will
  greatly simplify the calculations. This means we must pass in all 
  quantities we want to use from the calling routine which is fine.

*/
  time_dbg prj_drv_time = time_dbg("prj_drv") ;

  if( (opt/10) % 10 == 2){
     if ( (opt % 10) % 2 == 1){
       proj_HFB_real( com, intarr) ;
     } else {
       proj_HFB_cplx( com, intarr) ;
       }
  } else if( (opt/10) % 10 == 1){
    ;
//    proj_HF
    }

  prj_drv_time.end() ;

  return ;

} ;

void proj_HFB_real( common& com, std::vector<tei>& intarr){
/*
  Projected HFB for Real wavefunctions

    - This will set up all the necessary routines to do various projection.
    - Set up integration grids

    - Let's start with Number projection as defualt.  We will do repreated diagonalization 
      of the Number projected Hmailtonian

    - This is only Generalized right now

    - Matrices are complex

  H - Core hamiltonian
  oaoint - integrals in the orthogonal basis
  W - HFB wavefunction
  x - integration points
  wt - integration weights
*/

  int nbas = com.nbas() ;
  int maxit = com.mxscfit() ;
  int nele = com.nele() ;
  int nalp ;
  double mu = d0 ; 
  Eigen::MatrixXcd h ;
  Eigen::MatrixXcd xs, xsi ;
  Eigen::MatrixXcd R ;
  Eigen::MatrixXcd I, D ;
  Eigen::MatrixXcd p, k ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w ;
  time_dbg proj_HFB_real_time = time_dbg("proj_HFB_real") ;

  h.resize( 4*nbas, 4*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  D.resize( 2*nbas, 2*nbas) ;
  xs.resize( 2*nbas, 2*nbas) ;
  xsi.resize( 2*nbas, 2*nbas) ;
  p.resize( 2*nbas, 2*nbas) ;
  k.resize( 2*nbas, 2*nbas) ;
  w.moc.resize( 4*nbas, 4*nbas) ;
  w.eig.resize( 4*nbas) ;

  p.setZero() ;
  k.setZero() ;
  h.setZero() ;
  xs.setZero() ;
  xsi.setZero() ;
  xs.block( 0, 0, nbas, nbas).real() = com.getXS() ;
  xs.block( nbas, nbas, nbas, nbas) = xs.block( 0, 0, nbas, nbas) ;
  xsi = xs.inverse() ;
  h.block( 0, 0, nbas, nbas).real() = com.getH() ;
  h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
  transform( 0, xs, h.block( 0, 0, 2*nbas, 2*nbas)) ;
  h.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -h.block( 0, 0, 2*nbas, 2*nbas) ;

  /*
    Build rho and kappa from a previous HF calculation
  */

  nalp = nele/2 ;
  thermal_guess( nalp, nbas, p, k) ;
  w.moc.resize( 4*nbas, 4*nbas) ;
  w.eig.resize( 4*nbas) ;

/*
  Build the particle number integration grid
*/

  trapezoid* ngrid = new trapezoid( d0, d2*pi, 11) ;

  rrHFB_projection( nbas, h, intarr, w.moc, p, k, mu, ngrid, maxit, xs, xsi) ;

  proj_HFB_real_time.end() ;

  return ;

} ;

void proj_HFB_cplx( common& com, std::vector<tei>& intarr){
/*
  Projected HFB for Complex wavefunctions

    - This will set up all the necessary routines to do various projection.
    - Set up integration grids

    - Let's start with Number projection as defualt.  We will do repreated diagonalization 
      of the Number projected Hmailtonian

    - This is only Generalized right now

    - Matrices are complex

  H - Core hamiltonian
  oaoint - integrals in the orthogonal basis
  W - HFB wavefunction
  x - integration points
  wt - integration weights
*/

  int nbas = com.nbas() ;
  int nele = com.nele() ;
  int maxit = com.mxscfit() ;
  double mu = com.mu(), Ef = d0, djunk ; 
  cd energy ;
  Eigen::MatrixXd rtmp ;
  Eigen::MatrixXcd H ;
  Eigen::MatrixXcd m ;
  Eigen::MatrixXd s ;
  Eigen::MatrixXcd sc ;
  Eigen::MatrixXd Xs ;
  Eigen::MatrixXcd Xsc, R ;
  Eigen::MatrixXcd V ;
  Eigen::MatrixXcd U ;
  Eigen::MatrixXcd G ;
  Eigen::MatrixXcd D ;
  Eigen::MatrixXcd t ;
  Eigen::MatrixXcd rho, kappa ;
  Eigen::VectorXcd rl ( 2*nbas), ll ( 2*nbas) ;
  std::vector<tei> oaoint ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;
  time_dbg proj_HFB_cplx_time = time_dbg("proj_HFB_cplx") ;
//
  cd Uval, Vval ;
//

  rtmp.resize( nbas, nbas) ;
  H.resize( 2*nbas, 2*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  s.resize( 2*nbas, 2*nbas) ;
  sc.resize( 4*nbas, 4*nbas) ;
  Xs.resize( nbas, nbas) ;
  Xsc.resize( 2*nbas, 2*nbas) ;
  w.moc.resize( 2*nbas, 2*nbas) ;
  w.eig.resize( 2*nbas) ;

  m.resize( 4*nbas, 4*nbas) ; 
  R.resize( 2*nbas, 2*nbas) ;

  rho.resize( 2*nbas, 2*nbas) ;
  kappa.resize( 2*nbas, 2*nbas) ;
  V.resize( 2*nbas, 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;
  G.resize( nbas*2, nbas*2) ;
  D.resize( nbas*2, nbas*2) ;

/*
  Orthogonalize the 
    -one electron Hamiltonian
    -two electron integrals
    -the wavefunction.
*/  

  H.setZero() ;
  rtmp = com.getH() ;
//  print_mat( rtmp, " Hamiltonian matrix ") ;
  Xs = com.getXS() ;
//  print_mat( Xs, " Oao transform " ) ;
  H.block( 0, 0, nbas, nbas).real() =  Xs.adjoint()*rtmp*Xs ;
  H.block( nbas, nbas, nbas, nbas) = H.block( 0, 0, nbas, nbas) ;
/*
  Leave out chemical potential for now 
*/
  std::cout << " mu " << mu << std::endl ;
  H += mu*Eigen::MatrixXcd::Identity( 2*nbas, 2*nbas) ;
//  print_mat( H, " Oao H + mu ") ;

  oao( nbas, intarr, oaoint, Xs) ;

  s.setZero() ;
  s.block( 0, 0, nbas, nbas) = com.getS() ;
  s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
  print_mat( s, " s ") ;

  canort( s, Xsc, 2*nbas) ;

/*
  std::cout << " Checking the orthogonalization routine " << std::endl ;
  std::cout << " Xsc.adjoint()*s*Xsc " << std::endl ;
  std::cout << Xsc.adjoint()*s*Xsc << std::endl ;
*/

  load_wfn( w) ;

/*
  Orthogonalize the wavefunction 
*/
//  m.setZero() ;
//  m.block( 0, 0, 2*nbas, 2*nbas) = s ;
//  m.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = s ;

  U = w.moc ;
  w.moc = Xsc.adjoint()*s*U ;
//  w.moc.block( 0, 0, 2*nbas, 2*nbas) = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
//  U = w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
//  w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) = Xsc.adjoint()*s*U ;
//  w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas) = w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
/*
  std::cout << " Orthonormal " << std::endl ;
  std::cout << " W^{t}*W " << std::endl ;
*/
  R = w.moc.adjoint()*w.moc ;
  print_mat( R, " Orthogonal wfn") ;

/*
  Mix the homo and lumo
*/
  rl = ( w.moc.col( nele-1) + w.moc.col( nele))/std::sqrt(z2) ;
  ll = ( w.moc.col( nele-1) - w.moc.col( nele))/std::sqrt(z2) ;
  w.moc.col( nele-1) = rl ;
  w.moc.col( nele) = ll ;

  Ef = (w.eig( nele - 1) + w.eig( nele))/d2 ;
  for( int i=0; i < 2*nbas; i++){
    djunk = (w.eig( i) - Ef)/(kb*3.0e4) ;
    rl( i) = static_cast<cd>(d1/( d1 + std::exp( djunk))) ;
    }

  rho = U*U.adjoint() ;
  print_mat( rho, " 0 k rho") ;
  print_mat( rl, "occ #") ;
  U = w.moc  ;
  print_mat( U, "Mo cof") ;
  U = U*rl.asDiagonal() ;
  print_mat( U, "Thermalized") ;

  rho = w.moc*U.adjoint() ;

  print_mat( rho, " Thermalized density") ;
//  rho = w.moc.conjugate()*w.moc.transpose() ;

  kappa.block( 0, nbas, nbas, nbas) = rho.block( 0, 0, nbas, nbas) - rho.block( 0, 0, nbas, nbas)*rho.block( 0, 0, nbas, nbas) ;
  kappa.block( nbas, 0, nbas, nbas) = -kappa.block( 0, nbas, nbas, nbas) ;
//  kappa = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;
  std::cout << " Orthonormal rho and kappa " << std::endl ;
  print_mat( rho) ;
  print_mat( kappa) ;
/*  ctr2eg( oaoint, rho, G, nbas) ;
  ctrPairg( oaoint, kappa, D, nbas) ;
  t = (H + G/d2)*rho ;
  energy = t.trace() ;
  t = kappa.conjugate()*D ;
  energy += t.trace()/d4 ;
  std::cout << "    Electronic Energy: " << energy.real() << std::endl << std::endl ;
  std::cout << " rho before thermal " << std::endl ;
  print_mat( rho) ;

    Thermalize rho

  kappa.setZero() ;
  
  kappa.block( 0, nbas, nbas, nbas) = rho.block( 0, 0, nbas, nbas)*rho.block( 0, 0, nbas, nbas) - rho.block( 0, 0, nbas, nbas) ;
  kappa.block( nbas, 0, nbas, nbas) = -kappa.block( 0, nbas, nbas, nbas) ;
  print_mat( kappa) ;
  

  for( int qtz = 0; qtz < 2*nbas; qtz+=4){
    Vval = rl[qtz/4]*z1 ;
    Uval = (std::sqrt( z1 - std::pow( rl[qtz/4], 2)))*z1 ;
    U( qtz, qtz) = Uval ;
    U( qtz + 1, qtz + 1) = Uval ;
    V( qtz, qtz + 1) = Vval ;
    V( qtz + 1, qtz) = -Vval ;
    }

  for( int qtz = 2; qtz < 2*nbas; qtz+=4){
    Vval = rl[qtz/4 + nbas/2]*z1 ;
    Uval = (std::sqrt( z1 - std::pow( rl[qtz/4 + nbas/2], 2)))*z1 ;
    U( qtz, qtz) = Uval ;
    U( qtz + 1, qtz + 1) = Uval ;
    V( qtz, qtz + 1) = Vval ;
    V( qtz + 1, qtz) = -Vval ;
    }
*/

/*
  Build the particle number integration grid
*/

  trapezoid* ngrid = new trapezoid( d0, d2*pi, 11) ;

  cgHFB_projection_dbg( nbas, H, oaoint, w.moc, rho, kappa, ngrid, maxit) ; 

  proj_HFB_cplx_time.end() ;

  return ;

} ;

void rrHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, trapezoid*& ngrid, int& maxit, Eigen::Ref<Eigen::MatrixXcd> xs, Eigen::Ref<Eigen::MatrixXcd> xsi) {

/*
  Various debugging for the number projected HFB.
*/

  int iter = 0, iter_N = 0 ;
  int in ;
  int inmb = ngrid->ns() ;
  double b_ll, b_ul, nele = d4 ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  cd lshift = -z6, N = z0 ;
  cd z2pi ( d2*pi, d0), x ( 40.0e0, d0) ;
  cd vac_nrm, energy = z0 ; /* This is the vacuum normalization */
  cd s_theta, intE_g, intx_g, x_g, p_energy ;
  cd blah3, blah2, blah, fnx, e_theta, fnw;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::MatrixXcd t ;
  Eigen::MatrixXcd nX ;
  Eigen::MatrixXcd Rp ;
  Eigen::MatrixXcd Heff, H ;
  Eigen::MatrixXcd Sdotrho, Sdotkappa ;
  Eigen::MatrixXcd Hdotrho, Hdotkappa ;
  Eigen::MatrixXcd f_theta ;
  Eigen::MatrixXcd kappa_bar ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd Rrho ;
  Eigen::MatrixXcd Rkappa ;
  Eigen::MatrixXcd rho_i, kappa_i ;
  Eigen::MatrixXcd Rrho_i ;
  Eigen::MatrixXcd Rkappa_i ;
  Eigen::MatrixXcd tmp ;
  Eigen::MatrixXcd r_theta ;
  Eigen::MatrixXcd k_theta ;
  Eigen::MatrixXcd I, G, D, D_bar, Vv, lvlshft ;
  Eigen::MatrixXcd Zm1, M11, M22, mu ;
  Eigen::IOFormat tmpfmt(4, 0, ", ", "\n", "[", "]") ;

/*
  Eigen Vectors to accumulate quantities
*/

  Eigen::VectorXcd OHg( inmb), OSg( inmb), d_g( inmb) ;
  Eigen::MatrixXcd C_theta ;
  Eigen::MatrixXcd C_theta_i ;
//  std::ofstream F_OUT ;
  time_dbg cgHFB_projection_time = time_dbg("rrHFB_projection") ;
  t.resize( 2*nbas, 2*nbas) ;
  f_theta.resize( 2*nbas, 2*nbas) ;
  r_theta.resize( 2*nbas, 2*nbas) ;
  k_theta.resize( 2*nbas, 2*nbas) ;
  Zm1.resize( 2*nbas, 2*nbas) ;
  Vv.resize( 4*nbas, 4*nbas) ;
  tmp.resize( 2*nbas, 2*nbas) ;
/* Normal rho and kappa */
  p_rho.resize( 2*nbas, 2*nbas) ;
  p_kappa.resize( 2*nbas, 2*nbas) ;
/* Rotated rho and kappa */
  Rkappa.resize( 2*nbas, 2*nbas) ;
  Rrho.resize( 2*nbas, 2*nbas) ;
/* Inverse Rho and kappa  */
  kappa_i.resize( 2*nbas, 2*nbas) ;
  rho_i.resize( 2*nbas, 2*nbas) ;
/* Rotated Inverse  */
  Rkappa_i.resize( 2*nbas, 2*nbas) ;
  Rrho_i.resize( 2*nbas, 2*nbas) ;
//
  C_theta.resize( 2*nbas, 2*nbas) ;
  C_theta_i.resize( 2*nbas, 2*nbas) ;
  Sdotrho.resize( 2*nbas, 2*nbas) ;
  Sdotkappa.resize( 2*nbas, 2*nbas) ;
  Hdotrho.resize( 2*nbas, 2*nbas) ;
  Hdotkappa.resize( 2*nbas, 2*nbas) ;
  Rp.resize( 4*nbas, 4*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;
  H.resize( 4*nbas, 4*nbas) ;
  nX.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  I.setIdentity() ;
  Heff.setZero() ;
  Rp.setZero() ;
  G.resize( 2*nbas, 2*nbas) ;
  D.resize( 2*nbas, 2*nbas) ;
  D_bar.resize( 2*nbas, 2*nbas) ;
  M11.resize( 2*nbas, 2*nbas) ;
  M22.resize( 2*nbas, 2*nbas) ;
  lvlshft.resize( 4*nbas, 4*nbas) ;
  mu.resize( 4*nbas, 4*nbas) ;

  vac_nrm = z1 ;
  cd sgn = std::pow( -z1, 4*nbas*( 4*nbas + 1)/2) ;
  mu.setZero() ;
  mu.block( 0, 0, 2*nbas, 2*nbas) = -I ;
  mu.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I ;

  print_mat( rho, " initial rho") ;
  print_mat( kappa, " initial kappa") ;

  /*
    Put into the orthonormal basis and prepare level shifting density
  */

  transform( 1, xsi, rho) ;
  transform( 1, xsi, kappa) ;
  lvlshft.block( 0, 0, 2*nbas, 2*nbas) = rho ;
  lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
  lvlshft.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
  lvlshft.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
  p_rho = rho ;
  p_kappa = kappa ;

  do {
    /*
      Loop over the number projection grid and accumulate quantities 
    */
    rho_i = rho.inverse() ;
    kappa_i = kappa.inverse() ;

    /*
      Clear our accumulation matrices
    */

    Sdotrho.setZero() ;
    Sdotkappa.setZero() ;
    Hdotrho.setZero() ;
    Hdotkappa.setZero() ;

    intx_g = z0 ;
    intE_g = z0 ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      d_g(in) = ngrid->w()*std::exp( -zi*fnx*z4) ;

      /*
        Build our rotated quantities
      */

      nX = I*std::exp( zi*fnx) ;
      Rrho = nX*rho ;
      Rrho_i = Rrho.inverse() ;
      Rkappa = nX*kappa ;
      Rkappa_i = Rkappa.inverse() ;
      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
      C_theta = C_theta_i.inverse() ;
      r_theta = Rrho*C_theta*rho ;
      k_theta = Rrho*C_theta*kappa ;
      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;

      /* 
        Get the energy expression.
        Start by putting into the non orthogonal basis
      */

      transform( 1, xs, r_theta) ;
      transform( 1, xs, k_theta) ;
      transform( 1, xs, kappa_bar) ;

      ctr2eg( intarr, r_theta, G, nbas) ;
      ctrPairg( intarr, k_theta, D, nbas) ;
      ctrPairg( intarr, kappa_bar, D_bar, nbas) ;
      D /= z2 ;
      D_bar /= z2 ;
      transform( 1, xsi, r_theta) ;
      transform( 1, xsi, k_theta) ;
      transform( 1, xsi, kappa_bar) ;
      transform( 0, xs, G) ;
      transform( 0, xs, D) ;
      transform( 0, xs, D_bar) ;
      f_theta = h.block( 0, 0, 2*nbas, 2*nbas) + G ;
      t = (h.block( 0, 0, 2*nbas, 2*nbas) + f_theta)*r_theta ;
      OHg( in) = t.trace()/z2 ;
      t = kappa_bar*D/z2 ;
      OHg( in) += t.trace() ;

      /* 
        Get the overlap and the necessary matrix inverses 
      */

      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
      tmp = -Rrho*Rkappa_i.conjugate() + rho_i*kappa ;
      M11 = tmp.inverse() ;
      M22 = rho_i*kappa*(I - M11*rho_i*kappa) ;
      Zm1 = -Rkappa_i.conjugate()*M11*nX + M22*kappa_i ;
 
      if ( in == 0 ) {
        vac_nrm = std::abs(z1/pfaffian_H( Rp)) ;
        std::cout << " vac_nrm at iter = 0; in = 0 " << std::endl ;
        std::cout << vac_nrm << std::endl ;
        Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
        Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
        Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
        Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
        }

      std::cout << " sqrt(det(RP))" << std::endl ;
      blah = Rp.determinant() ;
      std::cout << vac_nrm*std::sqrt( blah) << std::endl ;
      blah = std::sqrt(nX.determinant()) ;
      blah2 = std::sqrt( rho.determinant()) ;
      blah3 = std::sqrt( C_theta_i.determinant()) ;
      energy = blah*blah3/blah2 ;
      std::cout << " Onishi Density Overlap" << std::endl ;
      std::cout << energy << std::endl ;

      OSg( in) = vac_nrm*sgn*pfaffian_H( Rp) ;
      std::cout << " pfaffian overlap  " << OSg( in) << std::endl ;

      /*
        Sdotrho
      */
      Sdotrho += d_g( in)*OSg( in)*Zm1/z2 ;

      /*
        Hdotrho
      */
      tmp = C_theta*rho*f_theta*nX ;
      tmp += f_theta*Rrho*C_theta ;
      tmp += -r_theta*f_theta*Rrho*C_theta ;
      tmp += -C_theta*rho*f_theta*r_theta*nX ;
/*  
*/  
      tmp += r_theta*D*Rkappa.conjugate()*C_theta/z2 ;
      tmp += -D*Rkappa.conjugate()*C_theta/z2 ;
      tmp += -C_theta*rho*D*kappa_bar*nX/z2 ;
/*  
*/  
      tmp += C_theta*kappa*D_bar*nX/z2 ;
      tmp += -k_theta*D_bar*Rrho*C_theta/z2 ;
      tmp += -C_theta*kappa*D_bar*r_theta*nX/z2 ;
      Hdotrho += d_g( in)*OSg( in)*(OHg( in)*Zm1/z2 + tmp) ;
/*  
  Sdotkappa
*/  
      Zm1 = Rkappa_i.conjugate()*M11*Rrho*kappa_i.conjugate() ;
      Sdotkappa += d_g( in)*OSg( in)*Zm1/z2 ;
/*  
  Hdotkappa
*/  
      tmp = C_theta*rho*f_theta*k_theta*nX.conjugate() ;
      tmp += C_theta*kappa*D_bar*k_theta*nX.conjugate()/z2 ;
      tmp += -C_theta*rho*D*nX.conjugate()/z2 ;
      tmp += -C_theta*rho*D*Rkappa.conjugate()*C_theta*kappa*nX.conjugate()/z2 ;
      Hdotkappa += d_g( in)*OSg( in)*(OHg( in)*Zm1 + tmp)/z2 ;
/*   
  Accumulate the overlap and projected energy
*/  
      intx_g += OSg( in)*d_g( in) ;
      intE_g += OSg( in)*OHg( in)*d_g( in) ;
      }

    ngrid->set_s() ;

    cnv_energy.push_back( intE_g/intx_g) ;

    tmp = Hdotrho/intx_g  - intE_g*Sdotrho/std::pow( intx_g, z2) ;
    Heff.block( 0, 0, 2*nbas, 2*nbas) = tmp ;
    tmp = Hdotkappa/intx_g - intE_g*Sdotkappa/std::pow( intx_g, z2) ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = tmp ;
    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;

    Heff += lvlshft*lshift ;

    /*
      Adjust the chemical potential until the density gives the proper number
      of particles.
    */

    iter_N = 0 ;
    b_ul = lambda ;
    b_ll = lambda ;

    /*
       Set some initial limits
    */

    do {
      b_ll -= d3 ;
      H = Heff + static_cast<cd>(b_ll)*mu ;
      H_diag.compute( H) ;
      Vv = H_diag.eigenvectors() ;
      rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
      N = rho.trace() ;
      } while ( static_cast<double>(N.real()) > static_cast<double>(nele)) ;

    std::cout << " Lower limit chemical potential " << b_ll << std::endl ;

    do {
      b_ul += d3 ;
      H = Heff + static_cast<cd>(b_ul)*mu ;
      H_diag.compute( H) ;
      Vv = H_diag.eigenvectors() ;
      rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
      N = rho.trace() ;
      } while ( static_cast<double>(N.real()) < static_cast<double>(nele)) ;

    std::cout << " upper limit chemical potential " << b_ul << std::endl ;

    H = Heff + static_cast<cd>(lambda)*mu ;

    while ( iter_N++ < 30) {
      H_diag.compute( H) ;
      Vv = H_diag.eigenvectors() ;
      rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
      N = rho.trace() ;
      if ( std::abs(static_cast<double>(N.real()) - static_cast<double>(nele)) < 1.0e-5){
        std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
        std::cout << "    chemical potential: " << lambda << std::endl ;
        std::cout << "    Particle Number : " << N.real() << std::endl << std::endl ;
        break ;
      } else {

/*
  Bisection method
*/

        if ( static_cast<double>(N.real()) - static_cast<double>(nele) < d0){
/*
  Too few electrons. Increase the chemical potential
*/

          b_ll = lambda ;
          lambda = (b_ul + b_ll)/d2 ;
        } else {
          b_ul = lambda ;
          lambda = (b_ul + b_ll)/d2 ;
          }
        }
      H = Heff + static_cast<cd>(lambda)*mu ;
      }

   kappa = Vv.block(  2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;
   std::cout << " rms difference in the density " << energy << std::endl ;
   if ( energy.real() < 1.0e-7 ) { break ;}

/*
   print_mat( rho, " rho ") ;
   print_mat( kappa, " kappa ") ;
*/
   lvlshft.block( 0, 0, 2*nbas, 2*nbas) = rho ;
   lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
   lvlshft.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
   lvlshft.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
   p_rho = rho ;
   p_kappa = kappa ;

//   tmp = Vv.block( 0, 2*nbas, 2*nbas, 2*nbas)*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint() ;
//   vac_nrm = std::sqrt( tmp.determinant()) ;
//   std::cout << " vac_nrm " << std::endl ;
//   std::cout << vac_nrm << std::endl ;

/*
  Check for convergence
*/  
//    intx_g = z0 ;
//    energy = z0 ;
//    kappa_i = kappa.inverse() ;
//    for ( in = 0; in < inmb; in++){
//      fnx = static_cast<cd>( ngrid->q()) ;
//      blah = ngrid->w()*std::exp( -zi*fnx*z4)/z2pi ;
///*
//  Build our rotated quantities
//*/
//      nX = I*std::exp( zi*fnx) ;
//      Rrho = nX*rho ;
//      Rkappa = nX*kappa ;
//      Rkappa_i = Rkappa.inverse() ;
//      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
//      C_theta = C_theta_i.inverse() ;
//      r_theta = Rrho*C_theta*rho ;
//      k_theta = Rrho*C_theta*kappa ;
//      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;
///* 
//  Get the energy expression
//*/
//      ctr2eg( intarr, r_theta, G, nbas) ;
//      ctrPairg( intarr, k_theta, D, nbas) ;
//      f_theta = h + G/z2 ;
//      t = f_theta*r_theta ;
//      blah2 = t.trace() ;
//      t = D*kappa_bar ;
//      blah2 += t.trace()/d4 ;
///* 
//  Get the overlap and the necessary matrix inverses 
//*/
//      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
//      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
//      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
//      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
//      blah3 = vac_nrm*sgn*pfaffian_H( Rp) ;
//      intx_g += blah*blah3 ;
//      energy += blah*blah2*blah3 ;
//    }

//    ngrid->set_s() ;

/*
  Energy
*/
//   energy = energy/intx_g ;
//   std::cout << "PHF Energy" << std::endl ;
//   std::cout << energy << std::endl ;
   
 } while ( ++iter <= 100) ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }
//
//  F_OUT.close() ;
//
  cgHFB_projection_time.end() ;

  return ;

} ;

void cgHFB_projection_dbg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, trapezoid*& ngrid, int& maxit) {

/*
  Various debugging for the number projected HFB.
*/

  int iter = 0 ;
  int in ;
  int inmb = ngrid->ns() ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  cd z2pi ( d2*pi, d0) ;
  cd vac_nrm, energy = z0 ; /* This is the vacuum normalization */
  cd blah3, blah2 = z0, blah = z0, e_theta, fnw, fnx, s_theta, intE_g = z0, intx_g = z0, x_g, p_energy=z0 ;
  std::vector<cd> cnv_den ;
  Eigen::MatrixXcd t ;
  Eigen::MatrixXcd nX ;
  Eigen::MatrixXcd Rp ;
  Eigen::MatrixXcd Heff ;
  Eigen::MatrixXcd Sdotrho ;
  Eigen::MatrixXcd Sdotkappa ;
  Eigen::MatrixXcd Hdotrho ;
  Eigen::MatrixXcd Hdotkappa ;
  Eigen::MatrixXcd f_theta ;
  Eigen::MatrixXcd kappa_bar ;
  Eigen::MatrixXcd p_rho ;
  Eigen::MatrixXcd p_kappa ;
  Eigen::MatrixXcd Rrho ;
  Eigen::MatrixXcd Rkappa ;
  Eigen::MatrixXcd rho_i ;
  Eigen::MatrixXcd kappa_i ;
  Eigen::MatrixXcd Rrho_i ;
  Eigen::MatrixXcd Rkappa_i ;
  Eigen::MatrixXcd tmp ;
  Eigen::MatrixXcd Rp_i ;
  Eigen::MatrixXcd r_theta ;
  Eigen::MatrixXcd k_theta ;
  Eigen::MatrixXcd I, G, D, D_bar, Vv ;
  Eigen::MatrixXcd Zm1, M11, M22, lvlshft ;

/*
  Eigen Vectors to accumulate quantities
*/

  Eigen::VectorXcd OHg( inmb), OSg( inmb), d_g( inmb) ;

/*
  C_theta = [-kappa*R^{*}kappa^{*} + rho*R*rho]^{-1}
*/

  Eigen::MatrixXcd C_theta ;
  Eigen::MatrixXcd C_theta_i ;
  time_dbg cgHFB_projection_time = time_dbg("cgHFB_projection") ;

  Rp_i.resize( 4*nbas, 4*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  f_theta.resize( 2*nbas, 2*nbas) ;
  r_theta.resize( 2*nbas, 2*nbas) ;
  k_theta.resize( 2*nbas, 2*nbas) ;
  Zm1.resize( 2*nbas, 2*nbas) ;
  Vv.resize( 4*nbas, 4*nbas) ;
  tmp.resize( 2*nbas, 2*nbas) ;
/* Normal rho and kappa */
  p_rho.resize( 2*nbas, 2*nbas) ;
  p_kappa.resize( 2*nbas, 2*nbas) ;
/* Rotated rho and kappa */
  Rkappa.resize( 2*nbas, 2*nbas) ;
  Rrho.resize( 2*nbas, 2*nbas) ;
/* Inverse Rho and kappa  */
  kappa_i.resize( 2*nbas, 2*nbas) ;
  rho_i.resize( 2*nbas, 2*nbas) ;
/* Rotated Inverse  */
  Rkappa_i.resize( 2*nbas, 2*nbas) ;
  Rrho_i.resize( 2*nbas, 2*nbas) ;

  C_theta.resize( 2*nbas, 2*nbas) ;
  C_theta_i.resize( 2*nbas, 2*nbas) ;
  Sdotrho.resize( 2*nbas, 2*nbas) ;
  Sdotkappa.resize( 2*nbas, 2*nbas) ;
  Hdotrho.resize( 2*nbas, 2*nbas) ;
  Hdotkappa.resize( 2*nbas, 2*nbas) ;
  Rp.resize( 4*nbas, 4*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;
  nX.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  I.setIdentity() ;
  Heff.setZero() ;
  Rp.setZero() ;
  G.resize( 2*nbas, 2*nbas) ;
  D.resize( 2*nbas, 2*nbas) ;
  D_bar.resize( 2*nbas, 2*nbas) ;
  M11.resize( 2*nbas, 2*nbas) ;
  M22.resize( 2*nbas, 2*nbas) ;
  lvlshft.resize( 4*nbas, 4*nbas) ;

/*
  I guess we will go on RMS density for convergence.
*/
  
//  print_mat( rho, "rho") ;
//  print_mat( kappa, " kappa ") ;

/*
  Get the proper sign for using the pfaffian.

  Vacuum normalization for quasi-particles is given by the expression sqrt( det(U^{t}U).
  Calculate it here for the initial iteration.

*/
  cd sgn = std::pow( -z1, 4*nbas*( 4*nbas + 1)/2) ;
//  std::cout << " Particle Number " << rho.trace() << std::endl ;
/*
  tmp = w.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint()*w.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
  vac_nrm = std::sqrt( tmp.determinant()) ;
  std::cout << " vac_nrm " << std::endl ;
  std::cout << vac_nrm << std::endl ;
  since we are populating from a HF determinant, we have no initial U to get the vacuum normalization
*/

  vac_nrm = z1 ; 

  do {
    /*
      Loop over the number projection grid and accumulate quantities
    */
    p_energy = energy ;
    p_rho = rho ;
    p_kappa = kappa ;
    rho_i = rho.inverse() ;
    kappa_i = kappa.inverse() ;
    /*
      Clear our accumulation matrices
    */
    Sdotrho.setZero() ;
    Sdotkappa.setZero() ;
    Hdotrho.setZero() ;
    Hdotkappa.setZero() ;
    intx_g = z0 ;
    intE_g = z0 ;
    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      d_g(in) = ngrid->w()*std::exp( -zi*fnx*z4)/z2pi ;

/*
  Build our rotated quantities
*/
      nX = I*std::exp( zi*fnx) ;
      Rrho = nX*rho ;
      Rrho_i = Rrho.inverse() ;
      Rkappa = nX*kappa ;
      Rkappa_i = Rkappa.inverse() ;
      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
      C_theta = C_theta_i.inverse() ;
      r_theta = Rrho*C_theta*rho ;
      k_theta = Rrho*C_theta*kappa ;
      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;
//
//      print_mat( nX, " Rotation matrix ") ;
//      print_mat( Rrho, " Rrho ") ;
//      print_mat( Rkappa, " Rkappa ") ;
//      print_mat( Rrho_i, " Rrho_i ") ;
//      print_mat( Rkappa_i, " Rkappa_i ") ;
//      print_mat( C_theta_i, " C_theta_i ") ;
//      std::cout << " det(C_theta_i) " << std::endl ;
//      std::cout << C_theta_i.determinant() << std::endl ;
//      print_mat( C_theta, " C_theta ") ; 
//      std::cout << " r_theta " << std::endl ;
//      print_mat( r_theta) ;
//      std::cout << " k_theta " << std::endl ;
//      print_mat( k_theta) ;
//      std::cout << " kappa_bar " << std::endl ;
//      print_mat( kappa_bar) ;
//
/* 
  Get the energy expression
*/
    ctr2eg( intarr, r_theta, G, nbas) ;
//    print_mat( G, " G") ;
    ctrPairg( intarr, k_theta, D, nbas) ;
//    print_mat( D, " D") ;
    ctrPairg( intarr, kappa_bar, D_bar, nbas) ;
//    print_mat( D_bar, " D_bar") ;
    f_theta = h + G/z2 ;
//    print_mat( f_theta. " f_theta") ;
    t = f_theta*r_theta ;
    OHg( in) = t.trace() ;
    t = D*kappa_bar ;
    OHg( in) += t.trace()/d4 ;
//    std::cout << " H " << std::endl ;
//    std::cout << OHg( in) << std::endl ;
    f_theta = h + G ;
/* 
  Get the overlap and the necessary matrix inverses 
*/
    Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
    Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
    Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
    Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
    tmp = -Rrho*Rkappa_i.conjugate() + rho_i*kappa ;
    M11 = tmp.inverse() ;
    M22 = rho_i*kappa*(I - M11*rho_i*kappa) ;
    Zm1 = -Rkappa_i.conjugate()*M11*nX + M22*kappa_i ;
    OSg( in) = vac_nrm*sgn*pfaffian_H( Rp) ;
/*
  Sdotrho
*/
    Sdotrho += d_g( in)*OSg( in)*Zm1/z2 ;
/*
  Hdotrho
*/
    tmp = C_theta*rho*f_theta*nX ;
    tmp += f_theta*Rrho*C_theta ;
    tmp += -r_theta*f_theta*Rrho*C_theta ;
    tmp += -C_theta*rho*f_theta*r_theta*nX ;
    tmp += r_theta*D*Rkappa.conjugate()*C_theta/z4 ;
    tmp += -D*Rkappa.conjugate()*C_theta/z4 ;
    tmp += -C_theta*rho*D*kappa_bar*nX/z4 ;
//
    tmp += C_theta*kappa*D_bar*nX/z4 ;
    tmp += -k_theta*D_bar*Rrho*C_theta/z4 ;
    tmp += -C_theta*kappa*D_bar*r_theta*nX/z4 ;
    Hdotrho += d_g( in)*OSg( in)*(OHg( in)*Zm1/z2 + tmp) ;
/*
  Sdotkappa
*/

    Zm1 = Rkappa_i.conjugate()*M11*Rrho*kappa_i.conjugate() ;
    Sdotkappa += d_g( in)*OSg( in)*(Zm1 - Zm1.transpose())/z4 ;
/*
  Hdotkappa
*/
    tmp = C_theta*rho*f_theta*k_theta*nX.conjugate() ;
    tmp += C_theta*kappa*D_bar*k_theta*nX.conjugate()/z4 ;
    tmp += -C_theta*rho*D*nX.conjugate()/z4 ;
    tmp += -C_theta*rho*D*Rkappa.conjugate()*C_theta*kappa*nX.conjugate()/z4 ;
    Hdotkappa += d_g( in)*OSg( in)*(OHg( in)*(Zm1 - Zm1.transpose())/z4 + (tmp - tmp.transpose())/z2) ;
/* 
  Accumulate the overlap and projected energy
*/
    intx_g += OSg( in)*d_g( in) ;
    intE_g += OSg( in)*OHg( in)*d_g( in) ;
    }

    ngrid->set_s() ;

    tmp = Hdotrho/intx_g  - intE_g*Sdotrho/std::pow( intx_g, z2) ;
    Heff.block( 0, 0, 2*nbas, 2*nbas) = (tmp + tmp.adjoint())/z2 ;
    tmp = Hdotkappa/intx_g - intE_g*Sdotkappa/std::pow( intx_g, z2) ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = (tmp - tmp.transpose())/z2 ;
    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
//    print_mat( Heff, " Heff") ;
/*
  Level shifting because how else will this work
*/

    lvlshft.block( 0, 0, 2*nbas, 2*nbas) = rho ;
    lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
    lvlshft.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
    lvlshft.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
    Heff += lvlshft/z4 ;

    H_diag.compute( Heff) ;
    Vv = H_diag.eigenvectors() ;
    print_mat( Vv, " eigenvectors") ;
    std::cout << "Are the vectors orthonormal" << std::endl ;
    std::cout << Vv.adjoint()*Vv << std::endl ;
    rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
    tmp = (rho + rho.adjoint())/z2 ;
    rho = tmp ;
//    print_mat( rho, " rho") ;
    kappa = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;
    tmp = (kappa - kappa.transpose())/z2 ;
    kappa = tmp ;
    tmp = Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
    vac_nrm = std::sqrt( tmp.determinant()) ;
//    std::cout << " vac_nrm " << std::endl ;
//    std::cout << vac_nrm << std::endl ;
//    print_mat( kappa, " kappa projected") ;
/*
  Check for convergence
*/  
    intx_g = z0 ;
    energy = z0 ;
    kappa_i = kappa.inverse() ;
    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      blah = ngrid->w()*std::exp( -zi*fnx*z4)/z2pi ;
/*
  Build our rotated quantities
*/
      nX = I*std::exp( zi*fnx) ;
      Rrho = nX*rho ;
      Rkappa = nX*kappa ;
      Rkappa_i = Rkappa.inverse() ;
      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
      C_theta = C_theta_i.inverse() ;
      r_theta = Rrho*C_theta*rho ;
      k_theta = Rrho*C_theta*kappa ;
      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;
/* 
  Get the energy expression
*/
      ctr2eg( intarr, r_theta, G, nbas) ;
      ctrPairg( intarr, k_theta, D, nbas) ;
      f_theta = h + G/z2 ;
      t = f_theta*r_theta ;
      blah2 = t.trace() ;
      t = D*kappa_bar ;
      blah2 += t.trace()/d4 ;
/* 
  Get the overlap and the necessary matrix inverses 
*/
      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
      blah3 = vac_nrm*sgn*pfaffian_H( Rp) ;
      intx_g += blah*blah3 ;
      energy += blah*blah2*blah3 ;
    }

    ngrid->set_s() ;

/*
  Energy
*/
   energy = energy/intx_g ;
   std::cout << "PHF Energy" << std::endl ;
   std::cout << energy << std::endl ;

   t = rho - p_rho ;
   energy = ( t*t.adjoint()).norm() ;
   t = kappa - p_kappa ;
   energy += ( t*t.adjoint()).norm() ;
   cnv_den.push_back( energy) ;
   if ( energy.real() < 0.00001 ) { break ;}
   } while ( ++iter <= 100) ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << std::endl ;
    }

  cgHFB_projection_time.end() ;

  return ;

} ;

cd pf_overlap( int& nbas, cd& c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, Eigen::Ref<Eigen::MatrixXcd> R){

/*
  nbas - number of basis functions
  c - normalization term
  rho - unrotated density
  kappa - unrotated abnormal density
  R - rotation matrix
*/

  cd pf ;
  Eigen::MatrixXcd kappa_i ;
  Eigen::MatrixXcd Rrho ;
  Eigen::MatrixXcd Rkappa ;
  Eigen::MatrixXcd Rkappa_i ;
  Eigen::MatrixXcd Rp ;
  Eigen::MatrixXcd I ;

  Rp.resize( 4*nbas, 4*nbas) ;
  kappa_i.resize( 2*nbas, 2*nbas) ;
  Rrho.resize( 2*nbas, 2*nbas) ;
  Rkappa.resize( 2*nbas, 2*nbas) ;
  Rkappa_i.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  I.setIdentity() ;
  Rrho = R*rho ;
  kappa_i = kappa.inverse() ;
  Rkappa = R*kappa ;
  Rkappa_i = Rkappa.inverse() ;

  Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
  Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
  Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
  Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;

  pf = pfaffian_A( Rp) ;

  Rkappa_i.resize( 0, 0) ;
  Rkappa.resize( 0, 0) ;
  Rrho.resize( 0, 0) ;
  kappa_i.resize( 0, 0) ;
  Rp.resize( 0, 0) ;

  return pf ;

} ;

