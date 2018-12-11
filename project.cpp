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
       proj_HFB_cplx( com , intarr) ;
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

  int i, j ;
  int nbas = com.nbas() ;
  int maxit = com.mxscfit() ;
  double djunk ,mu = com.mu() ;
  Eigen::MatrixXd rtmp ;
  Eigen::MatrixXcd Tmp ;
  Eigen::MatrixXcd H ;
  Eigen::MatrixXd s ;
  Eigen::MatrixXcd sc ;
  Eigen::MatrixXd Xs ;
  Eigen::MatrixXcd Xsc, R ;
  Eigen::MatrixXcd V ;
  Eigen::MatrixXcd U ;
  std::vector<tei> oaoint ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;
  time_dbg proj_HFB_cplx_time = time_dbg("proj_HFB_cplx") ;
  srand(time(NULL)) ;
//
  double tmp, dlt_occ ;
  std::vector<double> rl ;
  cd Uval, Vval ;
//

  rtmp.resize( nbas, nbas) ;
  Tmp.resize( 4*nbas, 4*nbas) ;
  H.resize( 2*nbas, 2*nbas) ;
  s.resize( 4*nbas, 4*nbas) ;
  sc.resize( 4*nbas, 4*nbas) ;
  Xs.resize( nbas, nbas) ;
  Xsc.resize( 4*nbas, 4*nbas) ;
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

  V.resize( 2*nbas, 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;

/*
  Orthogonalize the 
    -one electron Hamiltonian
    -two electron integrals
    -the wavefunction.
*/  

  H.setZero() ;
  rtmp = com.getH() ;
  Xs = com.getXS() ;
  H.block( 0, 0, nbas, nbas).real() =  Xs.adjoint()*rtmp*Xs ;
  H.block(  nbas, nbas, nbas, nbas) = H.block( 0, 0, nbas, nbas) ;

  oao( nbas, intarr, oaoint, Xs) ;

  s.setZero() ;
  s.block( 0, 0, nbas, nbas) = com.getS() ;
  s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
  s.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = s.block( 0, 0, 2*nbas, 2*nbas) ;

/*
  Generate a random number broken guess
    1. Build the diagonalized U and V.  We must have fractional
       occupation of each state otherwise rho will not be invertible.
    2. Use the eigenvalues from an HFB guess to populate the levels
*/

  load_wfn( w) ;

  U.setZero() ;
  V.setZero() ;

  for( int qtz=0; qtz < nbas; qtz++){
    djunk = (w.eig( qtz) - mu)/(kb*1.0e5) ;
    tmp = d1/( d1 + std::exp( djunk)) ;
    rl.push_back( tmp) ;
    }

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

  w.moc.block( 0, 0, 2*nbas, 2*nbas) = U ;
  w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas) = V ;
  w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) = V ;
  w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = U ;

/*
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





  canort( s, Xsc, 4*nbas ) ;

  Tmp = Xsc.adjoint()*s*w.moc ;
  w.moc = Tmp ;
 

  Build the particle number integration grid
*/

  trapezoid* ngrid = new trapezoid( d0, d2*pi, 11) ;

  cgHFB_projection( nbas, H, oaoint, w.moc, ngrid, maxit) ; 

  proj_HFB_cplx_time.end() ;

  return ;

} ;

void cgHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, trapezoid*& ngrid, int& maxit) {

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
  int in ;
  int inmb = ngrid->ns() ;
  double energy = d0 ;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> H_diag ;
  cd vac_nrm ; /* This is the vacuum normalization */
  cd blah2 = z0,blah = z0, e_theta, fnw, fnx, s_theta, intx_g = z0, x_g ;
  Eigen::MatrixXcd t ;
  Eigen::MatrixXcd nX ;
  Eigen::MatrixXcd Rp ;
  Eigen::MatrixXcd Heff ;
  Eigen::MatrixXcd Fdotrho ;
  Eigen::MatrixXcd Sdotrho ;
  Eigen::MatrixXcd SdotrhoE ;
  Eigen::MatrixXcd SdotrhoO ;
  Eigen::MatrixXcd Sdotkappa ;
  Eigen::MatrixXcd SdotkappaE ;
  Eigen::MatrixXcd SdotkappaO ;
  Eigen::MatrixXcd Hdotrho ;
  Eigen::MatrixXcd Hdotkappa ;
  Eigen::MatrixXcd f_theta ;
  Eigen::MatrixXcd rho ;
  Eigen::MatrixXcd kappa ;
  Eigen::MatrixXcd p_rho ;
  Eigen::MatrixXcd p_kappa ;
  Eigen::MatrixXcd Rrho ;
  Eigen::MatrixXcd Rkappa ;
  Eigen::MatrixXcd rho_i ;
  Eigen::MatrixXcd kappa_i ;
  Eigen::MatrixXcd Rrho_i ;
  Eigen::MatrixXcd Rkappa_i ;
  Eigen::MatrixXcd tmp ;
/*
  Rotated Densities
*/
  Eigen::MatrixXcd r_theta ;
  Eigen::MatrixXcd k_theta ;
  Eigen::MatrixXcd I, G, D, Vv ;
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
  kappa.resize( 2*nbas, 2*nbas) ;
  rho.resize( 2*nbas, 2*nbas) ;
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
  Fdotrho.resize( 2*nbas, 2*nbas) ;
  Sdotrho.resize( 2*nbas, 2*nbas) ;
  SdotrhoE.resize( 2*nbas, 2*nbas) ;
  SdotrhoO.resize( 2*nbas, 2*nbas) ;
  Sdotkappa.resize( 2*nbas, 2*nbas) ;
  SdotkappaE.resize( 2*nbas, 2*nbas) ;
  SdotkappaO.resize( 2*nbas, 2*nbas) ;
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

/*
  I guess we will go on RMS density for convergence.
*/
  
  rho = w.block( 2*nbas, 0, 2*nbas, 2*nbas).conjugate()*w.block( 2*nbas, 0, 2*nbas, 2*nbas).transpose() ;
  kappa = w.block( 2*nbas, 0, 2*nbas, 2*nbas).conjugate()*w.block( 0, 0, 2*nbas, 2*nbas).transpose() ;
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


  do {
    /*
      Loop over the number projection grid and accumulate quantities
    */
    p_rho = rho ;
    p_kappa = kappa ;
    tmp = w.block( 0, 0, 2*nbas, 2*nbas).adjoint()*w.block( 0, 0, 2*nbas, 2*nbas) ;
    vac_nrm = std::sqrt(tmp.determinant()) ;
    Sdotrho.setZero() ;
    SdotrhoE.setZero() ;
    SdotrhoO.setZero() ;
    Hdotrho.setZero() ;
    Hdotkappa.setZero() ;
    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      d_g(in) = ngrid->w()*std::exp( -zi*fnx*static_cast<cd>(d2)) ;
      /*
        Generate our Rotated quantities we will need.  This is memory hungry right now
        but one step at a time plz. The matrices are poorly conditioned so I will 
        use the Marquardt-Levenberg coefficient.  I will have to come up with a better
        system but for now lets debug.
      */
      nX = I*std::exp( zi*fnx) ;
      Rrho = nX*rho ;
      tmp = Rrho + I*1.0e-3 ;
      Rrho_i = tmp.inverse() ;
      Rkappa = nX*kappa ;
      tmp = Rkappa + I*1.0e-3 ; 
      Rkappa_i = tmp.inverse() ;
      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
      tmp = C_theta_i + I*1.0e-3 ;
      C_theta = tmp.inverse() ;
      r_theta = Rrho*C_theta*rho ;
      k_theta = Rrho*C_theta*kappa ;

/*
      std::cout << " rotation angle " << fnx << std::endl ;
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
*/
      /* Get the energy expression */
/*
      ctr2eg( intarr, r_theta, G, nbas) ;
      ctrPairg( intarr, k_theta, D, nbas) ;
      f_theta = h + G/d2 ;

      std::cout << " h " << std::endl ;
      print_mat( h) ;

      std::cout << " G " << std::endl ;
      print_mat( G) ;

      std::cout << " f_theta " << std::endl ;
      print_mat( f_theta) ;
*/
      t = f_theta*r_theta ;
      OHg( in) = t.trace() ;
      t = -k_theta.conjugate()*D ;
      OHg( in) += t.trace()/d4 ;

      /* Get the overlap */
      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = rho.conjugate()*kappa_i ;
      M22_i = Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).inverse() ;
      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
      M11_i = Rp.block( 0, 0, 2*nbas, 2*nbas).inverse() ;
      OSg( in) = vac_nrm*pfaffian( Rp) ;
      /* Scale the energy by the overlap */
      OHg( in) = OHg( in) ;
      x_g = d_g( in)*OSg( in) ;
/*      std::cout << " M22_i " << std::endl ;
      print_mat( M22_i) ;
      std::cout << " M11_i " << std::endl ;
      print_mat( M11_i) ; 
*/
//      /* Accumulate the overlap derivative for later*/
      Sdotrho = OSg( in)*( kappa_i*M22_i - Rkappa_i.conjugate()*M11_i*nX)/d2 ;
/*      std::cout << " Sdotrho " << std::endl ;
      print_mat( Sdotrho) ; 
*/
      SdotrhoO += d_g( in)*Sdotrho ;
      Sdotkappa =  OSg( in)*(kappa_i.conjugate()*nX.transpose()*M11_i*Rrho*kappa_i.conjugate() - kappa_i*M22_i*rho.conjugate()*kappa_i)/d2 ;
/*      std::cout << " Sdotkappa " << std::endl ;
      print_mat( Sdotkappa) ; 
*/
      SdotkappaO += d_g( in)*Sdotkappa ;
      /* 
        Get the hamiltonian derivatives accumulated wrt rho and kappa
      */
/*
      std::cout << " G " << std::endl ;
      print_mat( G) ;
*/
      f_theta += G/d2 ;
/*
      std::cout << " f_theta " << std::endl ;
      print_mat( f_theta) ;
*/
      tmp = C_theta*rho*f_theta*nX ;
      tmp += f_theta*Rrho*C_theta ;
      tmp += -Rrho*C_theta_i*rho*f_theta*Rrho*C_theta_i ;
      tmp += -C_theta_i*rho*f_theta*Rrho*C_theta_i*rho*nX ;
      Hdotrho += x_g*tmp ;
      tmp = x_g*C_theta*kappa*D.conjugate()*nX - x_g*Rrho*C_theta_i*kappa*D.conjugate()*Rrho*C_theta_i ;
      tmp += -C_theta_i*kappa*D.conjugate()*Rrho*C_theta_i*rho*nX ;
      Hdotrho -= x_g*(tmp + tmp.conjugate())/d4 ;
      /*
      */
      Hdotkappa += -x_g*Rkappa.conjugate()*C_theta_i*rho*f_theta*Rrho*C_theta_i ;
      Hdotkappa += -x_g*C_theta_i*rho*f_theta*Rrho*C_theta_i*kappa*nX.conjugate() ;
      tmp += -Rkappa.conjugate()*C_theta_i*kappa*D.conjugate()*Rrho*C_theta_i ;
      tmp += -C_theta_i*kappa*D.conjugate()*Rrho*C_theta_i*kappa*nX.conjugate() ;
      tmp += D.conjugate()*Rrho*C_theta ;
      Hdotkappa += x_g*(tmp + tmp.conjugate())/d8 ;

      /* 
        Accumulate the overlap and weight 
        int dg <0|R(g)|0>w(g)
        Scale the enregy by the overlap to genreate <0|H|g>
      */
      intx_g += OSg( in)*d_g( in) ;
      }

/*
    std::cout << "intx_g " << intx_g << std::endl ;
    std::cout << " SdotkappaO " << std::endl ;
    print_mat( SdotkappaO) ;
    std::cout << " SdotrhoO " << std::endl ;
    print_mat( SdotrhoO) ;
    std::cout << " Hdotrho " << std::endl ;
    print_mat(  Hdotrho) ;
    std::cout << " Hdotkappa " << std::endl ;
    print_mat(  Hdotkappa) ;
*/

    ngrid->set_s() ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      /*
        Generate the rotation matrix for particle number
      */
      nX = I*std::exp( zi*fnx) ;
      Rrho = nX*rho ;
      tmp = Rrho + I*1.0e-3 ;
      Rrho_i = tmp.inverse() ;
      Rkappa = nX*kappa ;
      tmp = Rkappa + I*1.0e-3 ; 
      Rkappa_i = tmp.inverse() ;
      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
      tmp = C_theta_i + I*1.0e-3 ;
      C_theta = tmp.inverse() ;
      r_theta = Rrho*C_theta*rho ;
      k_theta = Rrho*C_theta*kappa ;

      /* Get the overlap */
      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = rho.conjugate()*kappa_i ;
      M22_i = Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).inverse() ;
      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
      M11_i = Rp.block( 0, 0, 2*nbas, 2*nbas).inverse() ;
      /* Get the Overlap derivative accumulated */
      x_g = d_g( in)*OSg( in) ;
      Sdotrho = (kappa_i*M22_i - Rkappa_i*M11_i*nX)/d2 ;
      SdotrhoE += d_g( in)*(Sdotrho - OSg( in)*SdotrhoO/intx_g)*OHg( in) ;
      Sdotkappa = kappa_i.conjugate()*nX.transpose()*M11_i*Rrho*kappa_i.conjugate() - kappa_i*M22_i*rho.conjugate()*kappa_i ;
      SdotkappaE += d_g( in)*(Sdotkappa - OSg( in)*SdotkappaO/intx_g)*OHg( in) ;
      }

    ngrid->set_s() ;
    tmp = Hdotrho + SdotrhoE/intx_g ;
    Heff.block( 0, 0, 2*nbas, 2*nbas) = tmp ;
    tmp = Hdotkappa + SdotkappaE/intx_g ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = tmp ;;
/*  
    tmp = Hdotrho + SdotrhoE/intx_g ;
    Heff.block( 0, 0, 2*nbas, 2*nbas) = (tmp + tmp.adjoint())/d2 ;
    tmp = Hdotkappa + SdotkappaE/intx_g ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = (tmp - tmp.transpose())/d2 ; 
*/
    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;

    H_diag.compute( Heff) ;
    Vv = H_diag.eigenvectors().block( 0, 2*nbas, 4*nbas, 2*nbas) ;
    std::cout << " H_eff.eigenvectors()" << std::endl ;
    print_mat( Vv) ;
    rho = Vv.block( 2*nbas, 0, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 0, 2*nbas, 2*nbas).transpose() ;
    kappa = Vv.block( 2*nbas, 0, 2*nbas, 2*nbas).conjugate()*Vv.block( 0, 0, 2*nbas, 2*nbas).transpose() ;

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
    if ( std::real(energy) < 0.1 ) { break ;}
    } while ( ++iter <= maxit ) ;

  cgHFB_projection_time.end() ;

  return ;

} ;
