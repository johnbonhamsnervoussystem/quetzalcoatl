#include <cmath>
#include "common.h"
#include <complex>
#include "constants.h"
#include <Eigen/Dense>
#include "evalm.h"
#include <iostream>
#include "hfwfn.h"
#include "integr.h"
#include "project.h"
#include "qtzcntrl.h"
#include "solver.h"
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include <vector>
#include "wigner.h"
#include "wfn.h"

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

  int nbas = com.nbas() ;
  int maxit = com.mxscfit() ;
  Eigen::MatrixXd rtmp ;
  Eigen::MatrixXcd Tmp ;
  Eigen::MatrixXcd H ;
  Eigen::MatrixXd s ;
  Eigen::MatrixXcd sc ;
  Eigen::MatrixXd Xs ;
  Eigen::MatrixXcd Xsc ;
  std::vector<tei> oaoint ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;
  time_dbg proj_HFB_cplx_time = time_dbg("proj_HFB_cplx") ;

  rtmp.resize( nbas, nbas) ;
  Tmp.resize( 4*nbas, 2*nbas) ;
  H.resize( 2*nbas, 2*nbas) ;
  s.resize( 4*nbas, 4*nbas) ;
  sc.resize( 4*nbas, 4*nbas) ;
  Xs.resize( nbas, nbas) ;
  Xsc.resize( 4*nbas, 4*nbas) ;
  w.moc.resize( 4*nbas, 2*nbas) ;
  w.eig.resize( 4*nbas) ;

/*
  Eigen::MatrixXcd V ;
  Eigen::MatrixXcd R ;
  int nbas = com.nbas() ;
  Eigen::MatrixXcd m ;
  Eigen::MatrixXcd V ;
  Eigen::MatrixXcd U ;
  Eigen::MatrixXcd rho ;
  Eigen::MatrixXcd kappa ;
  Eigen::MatrixXcd Uinv ;

  m.resize( 4*nbas, 4*nbas) ;
  V.resize( 2*nbas, 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;

  V.resize( 2*nbas, 2*nbas) ;
  R.resize( 2*nbas, 2*nbas) ;
*/

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

  load_wfn( w) ;
  s.setZero() ;
  s.block( 0, 0, nbas, nbas) = com.getS() ;
  s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
  s.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = s.block( 0, 0, 2*nbas, 2*nbas) ;

  canort( s, Xsc, 4*nbas) ;
  sc.real() = s ;
  Tmp.block( 0, 0, 4*nbas, 2*nbas) = Xsc.adjoint()*sc*w.moc ;
  w.moc = Tmp.block( 0, 0, 4*nbas, 2*nbas) ;

/*
  Build the particle number integration grid
*/
  trapezoid* ngrid = new trapezoid( d0, d2*pi, 11) ;

/* 
   W should have
   
     V
     U
*/

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
      | V|
      | U|
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
  cd e_theta, fnw, fnx ;
  Eigen::MatrixXcd t ;
  Eigen::MatrixXcd nX ;
  Eigen::MatrixXcd Rp ;
  Eigen::MatrixXcd R ;
  Eigen::MatrixXcd rho ;
  Eigen::MatrixXcd kappa ;
  Eigen::MatrixXcd tmp ;
/*
  Rotated Densities
*/
  Eigen::MatrixXcd r_theta ;
  Eigen::MatrixXcd k_theta ;
  Eigen::MatrixXcd I ;
  Eigen::MatrixXcd G ;
  Eigen::MatrixXcd D ;
/*
  C_theta = [-kappa*R^{*}kappa^{*} + rho*R*rho]^{-1}
*/
  Eigen::MatrixXcd C_theta ;
  time_dbg cgHFB_projection_time = time_dbg("cgHFB_projection") ;

  tmp.resize( 2*nbas, 2*nbas) ;
  kappa.resize( 2*nbas, 2*nbas) ;
  rho.resize( 2*nbas, 2*nbas) ;
  C_theta.resize( 2*nbas, 2*nbas) ;
  Rp.resize( 4*nbas, 4*nbas) ;
  R.resize( 4*nbas, 4*nbas) ;
  nX.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  I.setIdentity() ;
  R.setZero() ;
  Rp.setZero() ;
  G.resize( 2*nbas, 2*nbas) ;
  D.resize( 2*nbas, 2*nbas) ;

  /*
    I guess we will go on RMS density for convergence.
  */
  
  rho = w.block( 0, 0 , 2*nbas, 2*nbas).conjugate()*w.block( 0, 0 , 2*nbas, 2*nbas).transpose() ;
  kappa = w.block( 0, 0, 2*nbas, 2*nbas).conjugate()*w.block( 2*nbas, 0 , 2*nbas, 2*nbas).transpose() ;

  do {
    /*
      Loop over the number projection grid
    */
    for ( in = 0; in < inmb; in++){
      std::cout << " in : " << in << std::endl ;
      fnx = static_cast<cd>( ngrid->q()) ;
      std::cout << " a " << std::endl ;
      fnw = ngrid->w()*std::exp( -zi*fnx) ;
      std::cout << " b " << std::endl ;
      /*
        Generate the rotation matrix for particle number
      */
      nX = I*std::exp( zi*fnx) ;
      /*
        Generate our Rotated quantities
      */
      tmp = rho*nX*rho - kappa*nX.conjugate()*kappa.conjugate() ;
      C_theta = tmp.inverse() ;
      r_theta = nX*rho*C_theta*rho ;
      k_theta = nX*rho*C_theta*kappa ;
      /*
        Get the energy for this state.
      */
      std::cout << " r_theta " << std::endl ;
      std::cout << r_theta << std::endl ;
      std::cout << " k_theta " << std::endl ;
      std::cout << k_theta << std::endl ;
      ctr2eg( intarr, r_theta, G, nbas) ;
      ctrPairg( intarr, k_theta, D, nbas) ;
      t = (h + G/d2)*r_theta ;
      std::cout << " (h + G)*r_theta " << std::endl ;
      std::cout << t.trace() << std::endl ;
      e_theta = t.trace() ;
      std::cout << " k " << std::endl ;
      t = k_theta.conjugate()*D ;
      std::cout << " l " << std::endl ;
      e_theta += t.trace()/d4 ;
      std::cout << " energy " << std::endl ;
      std::cout << e_theta << std::endl ;
    }

    ngrid->set_s() ;

    std::cout << " integral " << std::endl ;
    std::cout << energy << std::endl ;
    /*
      Check for convergence
    */
    t = Rp - R ;
    energy = (t*t.adjoint()).norm() ;
    std::cout << "  rms difference in the densities: " << energy << std::endl ;
    if ( std::real(energy) < d1) { break ;}
    } while ( ++iter <= maxit ) ;

  cgHFB_projection_time.end() ;

  return ;

} ;
