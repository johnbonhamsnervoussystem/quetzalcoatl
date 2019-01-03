#include <cmath>
#include "common.h"
#include <complex>
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "evalm.h"
#include <iostream>
#include "hfwfn.h"
#include "integr.h"
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

  if( true){
     if ( false){
//       proj_HFB_real() ;
     } else {
       proj_HFB_cplx( com, intarr) ;
       }
  } else {
    ;
//    proj_HF
    }

  prj_drv_time.end() ;

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
  int maxit = com.mxscfit() ;
  double djunk, mu = com.mu() ;
  cd energy ;
  Eigen::MatrixXd rtmp ;
  Eigen::MatrixXcd Tmp ;
  Eigen::MatrixXcd H ;
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
  std::vector<tei> oaoint ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;
  time_dbg proj_HFB_cplx_time = time_dbg("proj_HFB_cplx") ;
//
  double tmp, dlt_occ ;
  std::vector<double> rl ;
  cd Uval, Vval ;
//

  rtmp.resize( nbas, nbas) ;
  Tmp.resize( 2*nbas, 2*nbas) ;
  H.resize( 2*nbas, 2*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  s.resize( 2*nbas, 2*nbas) ;
  sc.resize( 4*nbas, 4*nbas) ;
  Xs.resize( nbas, nbas) ;
  Xsc.resize( 2*nbas, 2*nbas) ;
  w.moc.resize( 4*nbas, 4*nbas) ;
  w.eig.resize( 4*nbas) ;

/*
  Eigen::MatrixXcd V ;
  Eigen::MatrixXcd R ;
  int nbas = com.nbas() ;
  Eigen::MatrixXcd m ;*/
/*  Eigen::MatrixXcd rho ;
  Eigen::MatrixXcd kappa ;
  Eigen::MatrixXcd Uinv ;
  m.resize( 4*nbas, 4*nbas) ; */
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
  std::cout << " Hamiltonian matrix " << std::endl ;
  print_mat( rtmp) ;
  Xs = com.getXS() ;
  std::cout << " Oao transform " << std::endl ;
  print_mat( Xs) ;
  H.block( 0, 0, nbas, nbas).real() =  Xs.adjoint()*rtmp*Xs ;
  H.block( nbas, nbas, nbas, nbas) = H.block( 0, 0, nbas, nbas) ;
  H -= mu*Eigen::MatrixXcd::Identity( 2*nbas, 2*nbas) ;
  std::cout << " Oao H - mu " << std::endl ;
  print_mat( H) ;

  oao( nbas, intarr, oaoint, Xs) ;

  s.setZero() ;
  s.block( 0, 0, nbas, nbas) = com.getS() ;
  s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
  std::cout << " s " << std::endl ;
  print_mat( s) ;

  canort( s, Xsc, 2*nbas) ;

  std::cout << " Xsc.adjoint()*s*Xsc " << std::endl ;
  std::cout << Xsc.adjoint()*s*Xsc << std::endl ;

  load_wfn( w) ;

/*
  Orthogonalize the wavefunction 
*/

  U = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) ;
  w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = Xsc.adjoint()*s*U ;
  U = w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
  w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) = Xsc.adjoint()*s*U ;


  rho = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
  kappa = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;
  std::cout << " Orthonormal rho and kappa " << std::endl ;
  print_mat( rho) ;
  print_mat( kappa) ;
  ctr2eg( oaoint, rho, G, nbas) ;
  ctrPairg( oaoint, kappa, D, nbas) ;
  t = (H + G/d2)*rho ;
  energy = t.trace() ;
  t = kappa.conjugate()*D ;
  energy += t.trace()/d4 ;
  std::cout << "    Electronic Energy: " << energy.real() << std::endl << std::endl ;
  std::cout << " rho before thermal " << std::endl ;
  print_mat( rho) ;

/*
    Thermalize rho
*/

  kappa.setZero() ;
  for( int i=0; i < 2*nbas; i++){
    djunk = (w.eig( i) - mu)/(kb*1.0e6) ;
    tmp = d1/( d1 + std::exp( djunk)) ;
    rl.push_back( tmp) ;
    }

  for( int i=0; i < 2*nbas; i+=2){
    std::cout <<  U.row( i).dot( U.row( i)) << std::endl ;
    rho( i, i) = U.row( i).dot( U.row( i))*rl[i] ;
    rho( i+1, i+1) = U.row( i+1).dot( U.row( i+1))*rl[ i+1] ;
//    kappa( i, i+1) = std::sqrt(d1 - pow( rho(i,i), 2)) ;
//    kappa( i+1, i) = -kappa( i, i+1) ;
    }

  std::cout << " Thermalized density " << std::endl ;
  print_mat( rho) ;
  
  kappa.block( 0, nbas, nbas, nbas) = rho.block( 0, 0, nbas, nbas)*rho.block( 0, 0, nbas, nbas) - rho.block( 0, 0, nbas, nbas) ;
  kappa.block( nbas, 0, nbas, nbas) = -kappa.block( 0, nbas, nbas, nbas) ;
  print_mat( kappa) ;
  
/*
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
  w.moc.block( 0, 0, 2*nbas, 2*nbas) = U ;
  w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas) = V ;
  w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) = V ;
  w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = U ;

  rand_unitary( U) ;
  Xsc.setZero() ;
  Xsc.block( 0, 0, 2*nbas, 2*nbas) = U ;
  Xsc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = U.conjugate() ;
  Tmp = w.moc*Xsc ;
  w.moc = Tmp ;
  rand_unitary( U) ;
  Xsc.setZero() ;
  Xsc.block( 0, 0, 2*nbas, 2*nbas) = U ;
  Xsc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = U.conjugate() ;
  Tmp = Xsc*w.moc ;
  w.moc = Tmp ; 


  U = w.moc.block( 0, 0, 2*nbas, 2*nbas) ;
  V = w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas) ;
  std::cout << " U " << std::endl ;
  std::cout << U << std::endl ;
  std::cout << " V " << std::endl ;
  std::cout << V << std::endl ;

  std::cout << " U^{t}*U + V^{t}*V " << std::endl ;
  R = U.adjoint()*U + V.adjoint()*V ;
  std::cout << R << std::endl ;
  std::cout << " sqrt(det(U^{t}*U + V^{t}*V)) " << std::endl ;
  std::cout << std::sqrt(R.determinant()) << std::endl ;

  std::cout << " U^{T}*V + V^{T}*U " << std::endl ;
  std::cout << U.transpose()*V + V.transpose()*U << std::endl ;
  std::cout << " U^{T}*V + V^{T}*U " << std::endl ;
  std::cout << U*V.adjoint() + V.conjugate()*U.transpose() << std::endl ;


  Build the particle number integration grid
*/

  trapezoid* ngrid = new trapezoid( d0, d2*pi, 11) ;

  cgHFB_projection( nbas, H, oaoint, w.moc, rho, kappa, ngrid, maxit) ; 

  proj_HFB_cplx_time.end() ;


  return ;

} ;

void cgHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, trapezoid*& ngrid, int& maxit) {

/*
  Step one is to build a transition density.
    - BUild the rotation matrix for particle number
  
  input :
    nbas - number of basis functions
    h - core hamiltonian in the orthogonal basis
    intarr - vector of two electron integrals
    w - hfb wavefunction
      | U V^{*}|
      | V U^{*}|
    ngrid - a trapezoidal grid for integration
    maxit - the upper limit for the scf iterations

  local :
    iter - iterations through the scf loop
    t - local scratch
    Rp - Previous iteration density
    R - Current iteration density
*/

  int iter = 0 ;
  int in, zz = 0 ;
  int inmb = ngrid->ns() ;
  double energy = d0 ;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> H_diag ;
  cd vac_nrm ; /* This is the vacuum normalization */
  cd blah2 = z0,blah = z0, e_theta, fnw, fnx, s_theta, intE_g = z0, intx_g = z0, x_g ;
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
  Eigen::MatrixXcd Gd ;
/*
  Rotated Densities
*/
  Eigen::MatrixXcd r_theta ;
  Eigen::MatrixXcd k_theta ;
  Eigen::MatrixXcd I, G, D, D_bar, Vv ;
  Eigen::MatrixXcd M11_i, M22_i ;
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

  t.resize( 2*nbas, 2*nbas) ;
  f_theta.resize( 2*nbas, 2*nbas) ;
  r_theta.resize( 2*nbas, 2*nbas) ;
  k_theta.resize( 2*nbas, 2*nbas) ;
  Vv.resize( 4*nbas, 2*nbas) ;
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
  Gd.resize( 4*nbas, 4*nbas) ;

/*
  I guess we will go on RMS density for convergence.
*/
  
  rho_i = rho.inverse() ;
  kappa_i = kappa.inverse() ;

  std::cout << " rho " << std::endl ;
  print_mat( rho) ;

  std::cout << " kappa " << std::endl ;
  print_mat( kappa) ;

  std::cout << " Particle Number " << rho.trace() << std::endl ;

/*
  std::cout << " p*p - p " << std::endl ;
  std::cout << rho*rho - rho << std::endl ;

  std::cout << " -kappa*kappa^{t} " << std::endl ;
  std::cout << -kappa*kappa.adjoint() << std::endl ;

  std::cout << " p*p - kappa*kappa^{*} " << std::endl ;
  std::cout << -kappa*kappa.conjugate() + rho*rho  << std::endl ; 
*/

/*
  tmp = w.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint()*w.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
  print_mat( tmp) ;
  vac_nrm = std::sqrt(tmp.determinant()) ;
  std::cout << " vac_nrm " << std::endl ;
  std::cout << vac_nrm << std::endl ;
*/

  do {
    /*
      Loop over the number projection grid and accumulate quantities
    */
    p_rho = rho ;
    p_kappa = kappa ;
    Sdotrho.setZero() ;
    Sdotkappa.setZero() ;
    Hdotrho.setZero() ;
    Hdotkappa.setZero() ;
    cd sgn = std::pow( -z1, 4*nbas*( 4*nbas + 1)/2) ;
    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      std::cout << " rotation angle " << fnx << std::endl ;
      d_g(in) = ngrid->w()*std::exp( -zi*fnx*static_cast<cd>(d2)) ;
      std::cout << " weight " << d_g( in) << std::endl ;
      /*
        Generate our Rotated quantities we will need.  This is memory hungry right now
        but one step at a time plz. The matrices are poorly conditioned so I will 
        use the Marquardt-Levenberg coefficient.  I will have to come up with a better
        system but for now lets debug.
      */

      nX = I*std::exp( zi*fnx) ;
      std::cout << " Rotation matrix " << std::endl ;
      print_mat( nX) ;
      Rrho = nX*rho ;
      tmp = Rrho ;
      Rrho_i = tmp.inverse() ;
      Rkappa = nX*kappa ;
      tmp = Rkappa ;
      Rkappa_i = tmp.inverse() ;
      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
      tmp = C_theta_i ;
      C_theta = tmp.inverse() ;
      r_theta = Rrho*C_theta*rho ;
      k_theta = Rrho*C_theta*kappa ;
      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;


      std::cout << " Rrho " << std::endl ;
      print_mat( Rrho) ;
      std::cout << " Rkappa " << std::endl ;
      print_mat( Rkappa) ;
      std::cout << " Rrho_i " << std::endl ;
      print_mat( Rrho_i) ;
      std::cout << " Rkappa_i " << std::endl ;
      print_mat( Rkappa_i) ;
      std::cout << " C_theta_i " << std::endl ;
      print_mat( C_theta_i) ;
      std::cout << " det(C_theta_i) " << std::endl ;
      std::cout << C_theta_i.determinant() << std::endl ;
      std::cout << " C_theta " << std::endl ;
      print_mat( C_theta) ; 
      std::cout << " r_theta " << std::endl ;
      print_mat( r_theta) ;
      std::cout << " k_theta " << std::endl ;
      print_mat( k_theta) ;
      std::cout << " kappa_bar " << std::endl ;
      print_mat( kappa_bar) ;

      /* Get the energy expression */

      ctr2eg( intarr, r_theta, G, nbas) ;
      std::cout << " G " << std::endl ;
      print_mat( G) ;
      ctrPairg( intarr, k_theta, D, nbas) ;
      std::cout << " D " << std::endl ;
      print_mat( D) ;
      ctrPairg( intarr, kappa_bar, D_bar, nbas) ;
      std::cout << " D_bar " << std::endl ;
      print_mat( D_bar) ;
      f_theta = h + G/d2 ;
      std::cout << " h " << std::endl ;
      print_mat( h) ;

      std::cout << " f_theta " << std::endl ;
      print_mat( f_theta) ;

      t = f_theta*r_theta ;
      OHg( in) = t.trace() ;
      t = kappa_bar.transpose()*D ;
      OHg( in) = t.trace()/d4 ;
      std::cout << " Energy " << OHg( in) << std::endl ;
      f_theta = h + G ;

      /* Get the overlap */
      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = rho.conjugate()*kappa_i ;
//      M22_i = Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).inverse() ;
      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
//      M11_i = Rp.block( 0, 0, 2*nbas, 2*nbas).inverse() ;
      tmp = Rp.block( 0, 0, 2*nbas, 2*nbas) + Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).inverse() ;
      M11_i = tmp.inverse() ;
      tmp = Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) + Rp.block( 0, 0, 2*nbas, 2*nbas).inverse() ;
      M22_i = tmp.inverse() ;

      std::cout << " M22_i " << std::endl ;
      print_mat( M22_i) ;
      std::cout << " M11_i " << std::endl ;
      print_mat( M11_i) ; 

      OSg( in) = sgn*pfaffian( Rp) ;
      if ( in == 0){
        vac_nrm = OSg( in) ;
        }
      OSg( in) /= vac_nrm ;
      std::cout << " <O|g> " << std::endl ;
      std::cout << OSg( in) << std::endl ;
      /* Scale the energy by the overlap */
      OHg( in) = OHg( in) ;
      x_g = d_g( in)*OSg( in) ;

      /* Collect all the terms in h_{ij} */
      tmp = (-Rkappa_i.conjugate()*M11_i*nX + M22_i*kappa_i)/d2 ;
      Sdotrho += x_g*tmp ;
      tmp *= OHg( in) ;
      tmp += C_theta*rho*f_theta*nX ;
      tmp += f_theta*Rrho*C_theta ;
      tmp += -r_theta*f_theta*Rrho*C_theta ;
      tmp += -C_theta*rho*f_theta*r_theta*nX ;
      tmp += -D.transpose()*Rkappa.conjugate()*C_theta/d4 ;
      tmp += r_theta*D.transpose()*Rkappa.conjugate()*C_theta/d4 ;
      tmp += -C_theta*rho*D.transpose()*kappa_bar*nX/d4 ;

      tmp += C_theta*kappa*D_bar.transpose()*nX/d4 ;
      tmp += -k_theta*D_bar.transpose()*Rrho*C_theta/d4 ;
      tmp += -C_theta*kappa*D_bar.transpose()*r_theta*nX/d4 ;
      Hdotrho += tmp*x_g ;

      /* Collect all the terms in D_{ij} */
      tmp = Rkappa_i.conjugate()*M11_i*Rrho*Rkappa_i.conjugate()/d2 ;
      Sdotkappa += x_g*tmp ;
      tmp *= OHg( in) ;
      std::cout << " f_theta" << std::endl ;
      print_mat( f_theta) ;
      std::cout << " tmp 1" << std::endl ;
      print_mat( tmp) ;
      std::cout << " blah 1" << std::endl ;
      std::cout << kappa*nX.conjugate() << std::endl ;
      std::cout << " blah 2" << std::endl ;
      std::cout << C_theta*kappa*nX.conjugate() << std::endl ;
      std::cout << " blah 3" << std::endl ;
      std::cout << Rrho*C_theta*kappa*nX.conjugate() << std::endl ;
      std::cout << " blah 4" << std::endl ;
      std::cout << f_theta*Rrho*C_theta*kappa*nX.conjugate() << std::endl ;
      std::cout << " blah 5" << std::endl ;
      std::cout << rho*f_theta*Rrho*C_theta*kappa*nX.conjugate() << std::endl ;
      std::cout << " blah 6" << std::endl ;
      std::cout << C_theta*rho*f_theta*k_theta*nX.conjugate() << std::endl ;
      tmp += C_theta*rho*f_theta*k_theta*nX.conjugate() ;
      std::cout << " tmp 2" << std::endl ;
      print_mat( tmp) ;
      tmp += -C_theta*rho*D.transpose()*nX.conjugate()/d4 ;
      std::cout << " tmp 3" << std::endl ;
      print_mat( tmp) ;
      tmp += -C_theta*rho*D.transpose()*Rkappa.conjugate()*C_theta*kappa*nX.conjugate()/d4 ;
      std::cout << " tmp 4" << std::endl ;
      print_mat( tmp) ;
      tmp += C_theta*kappa*D_bar.transpose()*k_theta*nX.conjugate()/d4 ;
      std::cout << " tmp 5" << std::endl ;
      print_mat( tmp) ;
      Hdotkappa += tmp*x_g ;
      std::cout << " Hdotkappa  iteration" << zz++ << std::endl ;
      print_mat( Hdotkappa) ;

      /* Accumulate the overlap and projected energy*/
      intx_g += OSg( in)*d_g( in) ;
      intE_g += OSg( in)*OHg( in)*d_g( in) ;
      }

    ngrid->set_s() ;

    tmp = Hdotrho/intx_g - intE_g*Sdotrho/std::pow( intx_g, d2) ;
    std::cout << " Heff " << std::endl ;
    print_mat( tmp) ;
    Heff.block( 0, 0, 2*nbas, 2*nbas) = tmp ;
    tmp = Hdotkappa/intx_g - intE_g*Sdotkappa/std::pow( intx_g, d2) ;
    std::cout << " Deff " << std::endl ;
    print_mat( tmp) ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = tmp ;
    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
    Gd.block( 0, 0, 2*nbas, 2*nbas) = rho ;
    Gd.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
    Gd.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
    Gd.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;

    Heff += Gd/d2 ;

    H_diag.compute( Heff) ;
    Vv = H_diag.eigenvectors() ;
    std::cout << " H_eff.eigenvectors()" << std::endl ;
    print_mat( Vv) ;
    rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
    kappa = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;

    std::cout << " rho " << std::endl ;
    print_mat( rho) ;

    std::cout << " kappa " << std::endl ;
    print_mat( kappa) ;

    std::cout << " Particle Number " << rho.trace() << std::endl ;

    /*
      Check for convergence
    */

    t = rho - p_rho ;
    energy = (t*t.adjoint()).norm() ;
    t = kappa - p_kappa ;
    energy += (t*t.adjoint()).norm() ;
    std::cout << "  rms difference in the densities: " << energy << std::endl ;
//    if ( std::real(energy) < 0.00001 ) { break ;}
    if ( iter >= 0 ) { break ;}
    } while ( ++iter <= maxit ) ;

  cgHFB_projection_time.end() ;

  return ;

} ;
