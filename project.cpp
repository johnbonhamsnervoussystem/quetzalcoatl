#include <cmath>
#include "common.h"
#include <complex>
#include "constants.h"
#include <deque>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include "evalm.h"
#include <iostream>
#include <fstream>
#include "init_job.h"
#include "integr.h"
#include <iomanip>
#include "guess.h"
#include "hfrout.h"
#include <math.h>
#include "nbodyint.h"
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

void prj_drv( common& com){
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

*/

 /*
   Start with Ring and Shiekh implementation.  No logic needed.
 */
  proj_HFB_real( com) ;

  return ;

} ;

void proj_HFB_real( common& com){
/*
  Projected HFB for Real wavefunctions

    - This will set up all the necessary routines to do various projection.
    - Set up integration grids

    - Let's start with Number projection as defualt.  We will do repreated diagonalization 
      of the Number projected Hmailtonian

  H - Core hamiltonian
  W - HFB wavefunction
  x - integration points
  wt - integration weights
*/

  int nbas = com.nbas() ;
  int maxit = com.mxscfit() ;
  int nalp = com.nalp() ;
  int pngrid = com.ngrid() ;
  double Nalp = static_cast<double>( nalp) ;
  double nele = static_cast<double>( com.nele()) ;
  int Nele = com.nele() ;
  double nrm, mu = d0 ;
  cd norm, olap, cospi4, sinpi4 ;
  cd lshift = static_cast<cd>(com.lvlshft()) ;
  Eigen::MatrixXcd h ;
  Eigen::MatrixXd I ;
  Eigen::MatrixXcd p, k, D, U, V ;
  Eigen::MatrixXd u, v, a, b ;
  Eigen::VectorXcd homo, lumo ;
  Eigen::MatrixXcd t, mix, UV, scr ;
  nbodyint<Eigen::MatrixXcd>* X ;
  std::vector<tei>* r12int ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w ;
  time_dbg proj_HFB_real_time = time_dbg("proj_HFB_real") ;

/*
  scr - scratch space
*/
  t.resize( nbas, nbas) ;
  UV.resize( 2*nbas, nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  D.resize( 2*nbas, 2*nbas) ;
  p.resize( nbas, nbas) ;
  k.resize( nbas, nbas) ;
  u.resize( nbas, nbas) ;
  v.resize( nbas, nbas) ;
  mix.resize( 2, 2*nbas) ;
  homo.resize( 2*nbas) ;
  lumo.resize( 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;
  V.resize( 2*nbas, 2*nbas) ;
  scr.resize( 4*nbas, 4*nbas) ;
  a.resize( nbas, nbas) ;
  b.resize( nbas, nbas) ;
  w.moc.resize( 2*nbas, 2*nbas) ;
  w.eig.resize( 2*nbas) ;

  h.setZero() ;
  p.setZero() ;
  k.setZero() ;
  U.setZero() ;
  V.setZero() ;
  scr.setZero() ;
  I.setIdentity() ;
  a.setZero() ;
  b.setZero() ;

//  int blahhhh = 3 ;
//  cplx_SlaDet( com, blahhhh) ;
//  wavefunction_characterization( nbas, nele) ;
  
//  generate_pk( nbas, nalp, p, k) ;
  jorge_guess( p, k, Nalp, nrm) ;

/*
  Use the Hartree-Fock eigenvalues to build thermally occupied orbitals
*/

  trapezoid* ngrid = new trapezoid( d0, d2*pi, pngrid) ;

  initialize( 2, 1, com.hamil(), com, h, X, r12int, nbas) ;

  t = h.block( 0, 0, nbas, nbas) ;

  mu = d0 ;
  lshift = z4 ;

  ring_shiekh_rr( nbas, t, X, w.moc, w.eig, p, k, mu, lshift, ngrid, nele, maxit) ;

  print_mat( w.eig, " HFB Orbital Energies ") ;
  UV = w.moc.block( 0, nbas, 2*nbas, nbas) ;
  print_mat( UV, " HFB wavefunction ") ;

/*

  p.resize( 2*nbas, 2*nbas) ;
  k.resize( 2*nbas, 2*nbas) ;

  Fill the beta block of the HFB wavefunction

  V.block( nbas, nbas, nbas, nbas) = w.moc.block( 0, 0, nbas, nbas) ;
  U.block( nbas, 0, nbas, nbas) = w.moc.block( nbas, 0, nbas, nbas) ;

  Generate the UHF wavefunction by mixing the HOMO-LUMO in the alpha MOs

  homo = w.moc.col(nbas-1) ;
  lumo = w.moc.col(nbas) ;
  print_mat( homo, " homo ") ;
  print_mat( lumo, " lumo ") ;
  mix.row(0) = (homo + lumo)/z2 ;
  mix.row(1) = (homo - lumo)/z2 ;
  w.moc.col( nbas-1) = mix.row(0) ;
  w.moc.col( nbas) = mix.row(1) ;
  print_mat( w.moc, " mixed mos ") ;

  Now that the orbitals are mixed, fill the alpha component of the MOs.

  V.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nbas) ;
  U.block( 0, nbas, nbas, nbas) = -w.moc.block( nbas, 0, nbas, nbas) ;

  p = V.conjugate()*V.transpose() ;
  k = V.conjugate()*U.transpose() ;

  initialize( 2, 3, com.hamil(), com, h, X, r12int, nbas) ;

  We need to pass in the quasi particle normalization to properly calculate pfaffians.
  Using the unrotated density and pairing, calculate 1/pfaffian for the norm.

  w.moc.resize( 4*nbas, 4*nbas) ;
  w.eig.resize( 4*nbas) ;

  D = k.inverse() ;

  scr.block( 0, 0, 2*nbas, 2*nbas) = -p*D.conjugate() ;
  scr.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
  scr.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
  scr.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = D*p ;
  olap = pfaffian_H( scr) ;
  norm = z1/olap ;

  std::cout << " norm " << norm << std::endl ;

  mu = d0 ;
  lshift = z4 ;

  ring_shiekh_cg( nbas, h, X, w.moc, w.eig, p, k, mu, lshift, ngrid, nele, norm, maxit) ;

*/
/*
  Now that we have our number project UHF HFB wavefunction, use that as a starting point
  for the GHF initial guess using a rotation angle of pi/4

  cospi4 = std::cos( zpi/z4) ;
  sinpi4 = std::sin( zpi/z4) ;


  V.block( 0, 0, nbas, 2*nbas) = cospi4*w.moc.block( 0, 0, nbas, 2*nbas) + sinpi4*w.moc.block( nbas, 0, nbas, 2*nbas) ;
  V.block( nbas, 0, nbas, 2*nbas) = cospi4*w.moc.block( nbas, 0, nbas, 2*nbas) - sinpi4*w.moc.block( 0, 0, nbas, 2*nbas) ;

  U.block( 0, 0, nbas, 2*nbas) = cospi4*w.moc.block( 2*nbas, 0, nbas, 2*nbas) + sinpi4*w.moc.block( 3*nbas, 0, nbas, 2*nbas) ;
  U.block( nbas, 0, nbas, 2*nbas) = cospi4*w.moc.block( 3*nbas, 0, nbas, 2*nbas) - sinpi4*w.moc.block( 2*nbas, 0, nbas, 2*nbas) ;

  p = V.conjugate()*V.transpose() ;
  k = V.conjugate()*U.transpose() ;

  print_mat( p, " density ") ;
  print_mat( k, " pairing ") ;

  D.setRandom() ;
  UV = (D.adjoint() - D)/z2 ;

  ahm_exp( UV, D, 2*nbas, 0) ;

  D *= z1/z10 ;

  scr.block( 0, 0, 2*nbas, 2*nbas) = V ;
  scr.block( 2*nbas, 0, 2*nbas, 2*nbas) = U ;

  V = D*scr.block( 0, 0, 2*nbas, 2*nbas) ;
  U = D.conjugate()*scr.block( 2*nbas, 0, 2*nbas, 2*nbas) ;

  p = V.conjugate()*V.transpose() ;
  k = V.conjugate()*U.transpose() ;

  While the matrices have the correct properties, let's iron out any instabilities

  V = p ;
  p = (V.adjoint() + V)/z2 ;
  V = k ;
  k = (V.transpose() - V)/z2 ;

  print_mat( p, " density ") ;
  print_mat( k, " pairing ") ;

  D = k.inverse() ;

  scr.block( 0, 0, 2*nbas, 2*nbas) = -p*D.conjugate() ;
  scr.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
  scr.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
  scr.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = D*p ;
  olap = pfaffian_H( scr) ;
  norm = z1/olap ;

  std::cout << " norm " << norm << std::endl ;

  mu = d0 ;
  lshift = z4 ;

  ring_shiekh_cg( nbas, h, X, w.moc, w.eig, p, k, mu, lshift, ngrid, nele, norm, maxit) ;

*/

//  tr_ring_shiekh( nbas, h, X, V, U, mu, lshift, ngrid, nele, maxit) ;

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

//  int nbas = com.nbas() ;
//  int nele = com.nele() ;
//  int maxit = com.mxscfit() ;
//  double mu = com.mu(), Ef = d0, djunk ; 
//  cd energy ;
//  Eigen::MatrixXd rtmp ;
//  Eigen::MatrixXcd H ;
//  Eigen::MatrixXcd m ;
//  Eigen::MatrixXd s ;
//  Eigen::MatrixXcd sc ;
//  Eigen::MatrixXd Xs ;
//  Eigen::MatrixXcd Xsc, R ;
//  Eigen::MatrixXcd V ;
//  Eigen::MatrixXcd U ;
//  Eigen::MatrixXcd G ;
//  Eigen::MatrixXcd D ;
//  Eigen::MatrixXcd t ;
//  Eigen::MatrixXcd rho, kappa ;
//  Eigen::VectorXcd rl ( 2*nbas), ll ( 2*nbas) ;
//  std::vector<tei> oaoint ;
//  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;
//  time_dbg proj_HFB_cplx_time = time_dbg("proj_HFB_cplx") ;
//
//  cd Uval, Vval ;
//

//  rtmp.resize( nbas, nbas) ;
//  H.resize( 2*nbas, 2*nbas) ;
//  t.resize( 2*nbas, 2*nbas) ;
//  s.resize( 2*nbas, 2*nbas) ;
//  sc.resize( 4*nbas, 4*nbas) ;
//  Xs.resize( nbas, nbas) ;
//  Xsc.resize( 2*nbas, 2*nbas) ;
//  w.moc.resize( 2*nbas, 2*nbas) ;
//  w.eig.resize( 2*nbas) ;
//
//  m.resize( 4*nbas, 4*nbas) ; 
//  R.resize( 2*nbas, 2*nbas) ;
//
//  rho.resize( 2*nbas, 2*nbas) ;
//  kappa.resize( 2*nbas, 2*nbas) ;
//  V.resize( 2*nbas, 2*nbas) ;
//  U.resize( 2*nbas, 2*nbas) ;
//  G.resize( nbas*2, nbas*2) ;
//  D.resize( nbas*2, nbas*2) ;

/*
  Orthogonalize the 
    -one electron Hamiltonian
    -two electron integrals
    -the wavefunction.
*/  

//  H.setZero() ;
//  rtmp = com.getH() ;
//  print_mat( rtmp, " Hamiltonian matrix ") ;
//  Xs = com.getXS() ;
//  print_mat( Xs, " Oao transform " ) ;
//  H.block( 0, 0, nbas, nbas).real() =  Xs.adjoint()*rtmp*Xs ;
//  H.block( nbas, nbas, nbas, nbas) = H.block( 0, 0, nbas, nbas) ;
/*
  Leave out chemical potential for now 
*/
//  std::cout << " mu " << mu << std::endl ;
//  H += mu*Eigen::MatrixXcd::Identity( 2*nbas, 2*nbas) ;
//  print_mat( H, " Oao H + mu ") ;

//  oao( nbas, intarr, oaoint, Xs) ;
//
//  s.setZero() ;
//  s.block( 0, 0, nbas, nbas) = com.getS() ;
//  s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
//  print_mat( s, " s ") ;
//
//  canort( s, Xsc, 2*nbas) ;

/*
  std::cout << " Checking the orthogonalization routine " << std::endl ;
  std::cout << " Xsc.adjoint()*s*Xsc " << std::endl ;
  std::cout << Xsc.adjoint()*s*Xsc << std::endl ;
*/

//  load_wfn( w) ;

/*
  Orthogonalize the wavefunction 
*/
//  m.setZero() ;
//  m.block( 0, 0, 2*nbas, 2*nbas) = s ;
//  m.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = s ;

//  U = w.moc ;
//  w.moc = Xsc.adjoint()*s*U ;
//  w.moc.block( 0, 0, 2*nbas, 2*nbas) = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
//  U = w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
//  w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) = Xsc.adjoint()*s*U ;
//  w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas) = w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
/*
  std::cout << " Orthonormal " << std::endl ;
  std::cout << " W^{t}*W " << std::endl ;
*/
//  R = w.moc.adjoint()*w.moc ;
//  print_mat( R, " Orthogonal wfn") ;

/*
  Mix the homo and lumo
*/
//  rl = ( w.moc.col( nele-1) + w.moc.col( nele))/std::sqrt(z2) ;
//  ll = ( w.moc.col( nele-1) - w.moc.col( nele))/std::sqrt(z2) ;
//  w.moc.col( nele-1) = rl ;
//  w.moc.col( nele) = ll ;
//
//  Ef = (w.eig( nele - 1) + w.eig( nele))/d2 ;
//  for( int i=0; i < 2*nbas; i++){
//    djunk = (w.eig( i) - Ef)/(kb*3.0e4) ;
//    rl( i) = static_cast<cd>(d1/( d1 + std::exp( djunk))) ;
//    }
//
//  rho = U*U.adjoint() ;
//  print_mat( rho, " 0 k rho") ;
//  print_mat( rl, "occ #") ;
//  U = w.moc  ;
//  print_mat( U, "Mo cof") ;
//  U = U*rl.asDiagonal() ;
//  print_mat( U, "Thermalized") ;
//
//  rho = w.moc*U.adjoint() ;
//
//  print_mat( rho, " Thermalized density") ;
//  rho = w.moc.conjugate()*w.moc.transpose() ;

//  kappa.block( 0, nbas, nbas, nbas) = rho.block( 0, 0, nbas, nbas) - rho.block( 0, 0, nbas, nbas)*rho.block( 0, 0, nbas, nbas) ;
//  kappa.block( nbas, 0, nbas, nbas) = -kappa.block( 0, nbas, nbas, nbas) ;
//  kappa = w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//  std::cout << " Orthonormal rho and kappa " << std::endl ;
//  print_mat( rho) ;
//  print_mat( kappa) ;
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

  trapezoid* ngrid = new trapezoid( d0, d2*pi, 11) ;

  cgHFB_projection_dbg( nbas, H, oaoint, w.moc, rho, kappa, ngrid, maxit) ; 

  proj_HFB_cplx_time.end() ;

*/
  return ;

} ;

template < class matrix>
void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig_v, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, int& maxit) {

/*
  Various debugging for the number projected HFB.
*/

  int iter = 1, iter_N = 0 ;
  int in, inmb = ngrid->ns() ;
  double b_ll, b_ul ;
  double pnum = nele/d2 ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  cd N = z0, olap ;
  cd energy, pphase, nphase ;
  cd s_theta, intE_g, intx_g, fnx, w_phi, x_phi, y_phi ;
  cd etr, gtr, dtr, d_phi ;
  cd acc_etr, acc_gtr, acc_dtr ;
  cd shift ;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::MatrixXcd R, evec ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd DbaN, tmp, Y_tmp, I ;
  Eigen::MatrixXcd C_phi, G_phi, D_phi ;
  Eigen::MatrixXcd r_phi, k_phi, R_phi ;
  Eigen::MatrixXcd X_phi, Y_phi ;
  Eigen::MatrixXcd kbar_phi, t ;
  Eigen::MatrixXcd mu, mu_n, Heff, H ;

  time_dbg ring_shiekh_rr_projection_time = time_dbg("ring and shiekh rr") ;

/*
  Local Variables
    R - Generalized density matrix
    evec - Eigenvectors
    p_rho - density from previous iteration
    p_kappa - pairing from previous iteration
*/

  R.resize( 2*nbas, 2*nbas) ;
  evec.resize( 2*nbas, 2*nbas) ;
  p_rho.resize( nbas, nbas) ;
  p_kappa.resize( nbas, nbas) ;
  epsilonN.resize( nbas, nbas) ;
  gammaN.resize( nbas, nbas) ;
  lambdaN.resize( nbas, nbas) ;
  FockN.resize( nbas, nbas) ;
  DeltaN.resize( nbas, nbas) ;
  DbaN.resize( nbas, nbas) ;
  Y_tmp.resize( nbas, nbas) ;
  tmp.resize( nbas, nbas) ;
  I.resize( nbas, nbas) ;
  C_phi.resize( nbas, nbas) ;
  G_phi.resize( nbas, nbas) ;
  R_phi.resize( nbas, nbas) ;
  D_phi.resize( nbas, nbas) ;
  r_phi.resize( nbas, nbas) ;
  k_phi.resize( nbas, nbas) ;
  X_phi.resize( nbas, nbas) ;
  Y_phi.resize( nbas, nbas) ;
  kbar_phi.resize( nbas, nbas) ;
  t.resize( nbas, nbas) ;
  mu.resize( 2*nbas, 2*nbas) ;
  mu_n.resize( 2*nbas, 2*nbas) ;
  Heff.resize( 2*nbas, 2*nbas) ;
  H.resize( 2*nbas, 2*nbas) ;

  I.setIdentity() ;
  mu.setZero() ;
  mu.block( 0, 0, nbas, nbas) = -I ;
  mu.block( nbas, nbas, nbas, nbas) = I ;
  mu_n.setIdentity() ;

  /*
    Prepare the generalized density with the initial guess
  */
  R.block( 0, 0, nbas, nbas) = rho ;
  R.block( 0, nbas, nbas, nbas) = kappa ;
  R.block( nbas, 0, nbas, nbas) = kappa.adjoint() ;
  R.block( nbas, nbas, nbas, nbas) = I - rho.conjugate() ;

  do {
    /*
      Save the previous matrices
    */
    p_rho = rho ;
    p_kappa = kappa ;

    /*
      Clear our accumulation matrices for the Hamiltonian
    */

    FockN.setZero() ;
    DeltaN.setZero() ;
    Y_tmp.setZero() ;

    intx_g = z0 ;
    intE_g = z0 ;
    acc_etr = z0 ;
    acc_gtr = z0 ;
    acc_dtr = z0 ;

/*
  Generate int d(phi) <O|phi> dphi and integral of the overlap derivative.
*/

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      w_phi = ngrid->w() ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>(nele))/( z2*zpi) ;
      pphase = std::exp( z2*zi*fnx) ;

      tmp = I + rho*( pphase - z1) ;
      C_phi = pphase*tmp.inverse() ;
      olap =  std::exp( zi*z2*fnx*static_cast<cd>(nbas))/C_phi.determinant() ;

      x_phi = d_phi*olap ;
      Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

      intx_g += w_phi*x_phi ;

      }

    ngrid->set_s() ;

    Y_tmp /= intx_g ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      w_phi = ngrid->w() ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>(nele))/( z2*zpi) ;
      pphase = std::exp( z2*zi*fnx) ;
      nphase = std::exp( -z2*zi*fnx) ;

/*
  Build our rotated quantities
*/

      tmp = I + rho*( pphase - z1) ;
      C_phi = pphase*tmp.inverse() ;
      r_phi = C_phi*rho ;
      k_phi = C_phi*kappa ;
      kbar_phi = pphase*kappa*C_phi.conjugate() ;
      olap =  std::exp( zi*z2*fnx*static_cast<cd>(nbas))/C_phi.determinant() ;
      x_phi = d_phi*olap ;
      y_phi = x_phi/intx_g ;
      Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi - Y_tmp ;

/*
  Calculate the Coulomb and Pairing Matrices
*/

      W->contract( r_phi, k_phi) ;
      evec = W->getG() ;
      G_phi = evec.block( 0, 0, nbas, nbas) ;
      D_phi = evec.block( 0, nbas, nbas, nbas) ;

/*
  Generate the one and two body hamiltonian elements
*/

      t = h*r_phi ;
      etr = z2*t.trace() ;
      epsilonN = w_phi*y_phi*(Y_phi*etr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi) ;
      t = G_phi*r_phi ;
      gtr = t.trace() ;
      gammaN = w_phi*y_phi*(Y_phi*gtr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*G_phi*C_phi) ;
      t = D_phi*kbar_phi.adjoint() ;
      dtr = -t.trace()/z2 ;
      t = D_phi.transpose()*kbar_phi.conjugate() ;
      dtr -= t.trace()/z2 ;
      lambdaN = -w_phi*y_phi*(Y_phi*dtr + z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint()) ;

      FockN += epsilonN + gammaN + lambdaN ;
      DeltaN += w_phi*nphase*y_phi*C_phi*D_phi ;

/*
  Accumulate the projected energy
*/
      acc_etr += w_phi*y_phi*etr ;
      acc_gtr += w_phi*y_phi*gtr ;
      acc_dtr -= w_phi*y_phi*dtr ;
      intE_g += w_phi*y_phi*( etr + gtr - dtr) ;
      }

    std::cout << " acc_etr " << acc_etr << std::endl ;
    std::cout << " acc_gtr " << acc_gtr << std::endl ;
    std::cout << " acc_dtr " << acc_dtr << std::endl ;

    ngrid->set_s() ;

    cnv_energy.push_back( intE_g) ;

/*
  Symmetrize the normal density component
*/

    t = (FockN + FockN.adjoint())/z2 ;
    Heff.block( 0, 0, nbas, nbas) = t ;

/*
  Symmetrize the pairing component.  Since kab = -kba we do
  DN + DN^{T} rather than DN - DN^{T}
*/

    t = (DeltaN + DeltaN.transpose())/z2 ;
    Heff.block( 0, nbas, nbas, nbas) = t ;

    Heff.block( nbas, nbas, nbas, nbas) = -Heff.block( 0, 0, nbas, nbas).conjugate() ;
    Heff.block( nbas, 0, nbas, nbas) = Heff.block( 0, nbas, nbas, nbas).adjoint() ;

    if ( false ){
//    if ( iter < 5 ){
/*
  Slowly scale down the level shifting
*/
      shift = lshift/static_cast<cd>(iter) ;
      }

//    Heff += -shift*R ;

  /*
    Adjust the chemical potential until the density gives the proper number
    of particles.
  */

    if ( false ){
      iter_N = 0 ;

      b_ul = lambda ;
      b_ll = lambda ;

/*
  Do the bisection method to find the value of the chemical potential.
*/

      do {
        b_ll -= d3 ;
        H = Heff + static_cast<cd>(b_ll)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 2*nbas, nbas)*evec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        } while ( static_cast<double>(N.real()) > pnum) ;

      std::cout << " Lower limit chemical potential " << b_ll << std::endl ;
      iter_N = 0 ;

      do {
        b_ul += d3 ;
        H = Heff + static_cast<cd>(b_ul)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 2*nbas, nbas)*evec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        } while ( static_cast<double>(N.real()) < pnum) ;

      std::cout << " upper limit chemical potential " << b_ul << std::endl ;
      iter_N = 0 ;

      lambda = ( b_ll + b_ul)/d2 ;

      while ( iter_N++ < 30) {
        H = Heff + static_cast<cd>(lambda)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 2*nbas, nbas)*evec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        if ( std::abs(static_cast<double>(N.real()) - pnum) < 1.0e-5){
          std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
          std::cout << "    chemical potential: " << lambda << std::endl ;
          std::cout << "    Particle Number : " << N.real() << std::endl << std::endl ;
          break ;
        } else {

          if ( static_cast<double>(N.real()) - pnum < d0){
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
      }
    } else {
/*
  If there is no chemical potential then simply diagonalize the Hamiltonian
*/
      H_diag.compute( Heff) ;
      evec = H_diag.eigenvectors() ;
      R = evec.block( 0, 0, 2*nbas, nbas)*evec.block( 0, 0, 2*nbas, nbas).adjoint() ;
      rho = R.block( 0, 0, nbas, nbas) ;
      }
/*
  Grab kappa
*/

   kappa = R.block( 0, nbas, nbas, nbas) ;

/*
  Check for convergence
*/

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;

   std::cout << " rms difference in the density " << energy << std::endl ;
   if ( energy.real() < 1.0e-7) {
     c = H_diag.eigenvectors() ;
     R.col(0) = H_diag.eigenvalues() ;
     eig_v = R.col(0).real() ;
     std::cout << " converged in the density " << std::endl ;
     break ;
     }
 
 } while ( ++iter < 250) ;

 for( unsigned int xq=0; xq < cnv_den.size(); xq++){
   std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
   }

  ring_shiekh_rr_projection_time.end() ;

  return ;

} ;

template void ring_shiekh_rr( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::VectorXd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd, trapezoid*&, double&, int&) ;

template < class matrix>
void ring_shiekh_cg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, cd &norm, int& maxit) {

/*
  Complex Generalized Number projected HFB
*/

  int iter = 1, iter_N = 0 ;
  int in, inmb = ngrid->ns() ;
  double b_ll, b_ul ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  cd N = z0, olap ;
  cd nrm = norm ;
  cd energy, pphase, nphase, shift ;
  cd s_theta, intE_g, intx_g, fnx, w_phi, x_phi, y_phi ;
  cd etr, gtr, dtr, d_phi ;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::VectorXcd teig ;
  Eigen::MatrixXcd R, evec ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd tmp, Y_tmp, I ;
  Eigen::MatrixXcd C_phi, G_phi, D_phi ;
  Eigen::MatrixXcd r_phi, k_phi, R_phi ;
  Eigen::MatrixXcd X_phi, Y_phi ;
  Eigen::MatrixXcd kbar_phi, t, t1 ;
  Eigen::MatrixXcd mu, mu_n, Heff, H ;

  time_dbg ring_shiekh_cg_time = time_dbg("ring and shiekh cg") ;

/*
  Local :
    R - Generalized density
    evec - eigenvectors/scratch space
*/

  R.resize( 4*nbas, 4*nbas) ;
  eig.resize( 4*nbas) ;
  evec.resize( 4*nbas, 4*nbas) ;
  p_rho.resize( 2*nbas, 2*nbas) ;
  p_kappa.resize( 2*nbas, 2*nbas) ;
  epsilonN.resize( 2*nbas, 2*nbas) ;
  gammaN.resize( 2*nbas, 2*nbas) ;
  lambdaN.resize( 2*nbas, 2*nbas) ;
  FockN.resize( 2*nbas, 2*nbas) ;
  DeltaN.resize( 2*nbas, 2*nbas) ;
  Y_tmp.resize( 2*nbas, 2*nbas) ;
  tmp.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  C_phi.resize( 2*nbas, 2*nbas) ;
  G_phi.resize( 2*nbas, 2*nbas) ;
  R_phi.resize( 2*nbas, 2*nbas) ;
  D_phi.resize( 2*nbas, 2*nbas) ;
  r_phi.resize( 2*nbas, 2*nbas) ;
  k_phi.resize( 2*nbas, 2*nbas) ;
  X_phi.resize( 2*nbas, 2*nbas) ;
  Y_phi.resize( 2*nbas, 2*nbas) ;
  kbar_phi.resize( 2*nbas, 2*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  t1.resize( 2*nbas, 2*nbas) ;
  mu.resize( 4*nbas, 4*nbas) ;
  mu_n.resize( 4*nbas, 4*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;
  H.resize( 4*nbas, 4*nbas) ;

  I.setIdentity() ;
  mu.setZero() ;
  mu.block( 0, 0, 2*nbas, 2*nbas) = -I ;
  mu.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I ;
  mu_n.setIdentity() ;

  /*
    Prepare level shifting density
  */

  R.block( 0, 0, 2*nbas, 2*nbas) = rho ;
  R.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
  R.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
  R.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;

  do {
    /*
      Loop over the number projection grid and accumulate quantities 
      Clear our accumulation matrices
    */

    p_rho = rho ;
    p_kappa = kappa ;
 
    FockN.setZero() ;
    DeltaN.setZero() ;
    Y_tmp.setZero() ;

    intx_g = z0 ;
    intE_g = z0 ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      w_phi = ngrid->w() ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
      pphase = std::exp( z2*zi*fnx) ;

      R_phi = I*std::exp( zi*fnx) ;
      tmp = I + rho*( pphase - z1) ;
      C_phi = pphase*tmp.inverse() ;

      t = R_phi*kappa ;
      t1 = t.inverse() ;
      evec.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*t1.conjugate() ;
      evec.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      evec.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      evec.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa.inverse()*rho ;
      olap = nrm*pfaffian_H( evec) ;

      x_phi = d_phi*olap ;
      Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

      intx_g += w_phi*x_phi ;

      }

    ngrid->set_s() ;

    Y_tmp /= intx_g ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      w_phi = ngrid->w() ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
      pphase = std::exp( z2*zi*fnx) ;
      nphase = std::exp( -z2*zi*fnx) ;

/*
  Build our rotated quantities
*/

      tmp = I + rho*( pphase - z1) ;
      C_phi = pphase*tmp.inverse() ;
      r_phi = C_phi*rho ;
      k_phi = C_phi*kappa ;
      kbar_phi = pphase*kappa*C_phi.conjugate() ;
      R_phi = I*std::exp( zi*fnx) ;
      t = R_phi*kappa ;
      t1 = t.inverse() ;
      evec.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*t1.conjugate() ;
      evec.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      evec.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      evec.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa.inverse()*rho ;
      olap = nrm*pfaffian_H( evec) ;
      x_phi = d_phi*olap ;
      y_phi = x_phi/intx_g ;
      Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi - Y_tmp ;

      W->contract( r_phi, k_phi) ;
      evec = W->getG() ;
      G_phi = evec.block( 0, 0, 2*nbas, 2*nbas) ;
      D_phi = evec.block( 0, 2*nbas, 2*nbas, 2*nbas) ;

      t = h*r_phi ;
      etr = t.trace() ;
      epsilonN = w_phi*y_phi*(Y_phi*etr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi) ;
      t = G_phi*r_phi ;
      gtr = t.trace() ;
      gammaN = w_phi*y_phi*(Y_phi*gtr/z2 + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*G_phi*C_phi) ;
      t = D_phi*kbar_phi.adjoint() ;
      dtr = -t.trace()/z2 ;
      t = D_phi.transpose()*kbar_phi.conjugate() ;
      dtr -= t.trace()/z2 ;
      lambdaN = -w_phi*y_phi*(Y_phi*dtr/z2 + z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint()) ;

      FockN += epsilonN + gammaN + lambdaN ;
      DeltaN += w_phi*nphase*y_phi*C_phi*D_phi ;

/*
  Accumulate the overlap and projected energy
*/
      intE_g += w_phi*y_phi*( etr + gtr/z2 - dtr/z2) ;
      }

    ngrid->set_s() ;

    cnv_energy.push_back( intE_g) ;

    t = (FockN + FockN.adjoint())/z2 ;
    FockN = t ;
    t = (DeltaN - DeltaN.transpose())/z2 ;
    DeltaN = t ;

    Heff.block( 0, 0, 2*nbas, 2*nbas) = FockN ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = DeltaN ;
    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint() ;

    if ( iter < 5 ){
/*
  Slowly scale down the level shifting
*/
      shift = lshift/static_cast<cd>(iter) ;
      }

    Heff += -shift*R ;
/*
  Adjust the chemical potential until the density gives the proper number
  of particles.
*/

    if ( true ) {
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
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
        rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
        N = rho.trace() ;
        } while ( static_cast<double>(N.real()) > nele) ;
    
      std::cout << " Lower limit chemical potential " << b_ll << std::endl ;
      iter_N = 0 ;
    
      do {
        b_ul += d3 ;
        H = Heff + static_cast<cd>(b_ul)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
        rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
        N = rho.trace() ;
        } while ( static_cast<double>(N.real()) < nele) ;
    
      std::cout << " upper limit chemical potential " << b_ul << std::endl ;
      iter_N = 0 ;
    
      lambda = ( b_ll + b_ul)/d2 ;
    
      while ( iter_N++ < 30) {
        H = Heff + static_cast<cd>(lambda)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
        rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
        N = rho.trace() ;
        if ( std::abs(static_cast<double>(N.real()) - nele) < 1.0e-5){
          std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
          std::cout << "    chemical potential: " << lambda << std::endl ;
          std::cout << "    Particle Number : " << N.real() << std::endl << std::endl ;
          break ;
        } else {
  
          if ( static_cast<double>(N.real()) - nele < d0){
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
        }
   } else {
/*
  If there is no chemical potential then simply diagonalize the Hamiltonian
*/
     H_diag.compute( Heff) ;
     evec = H_diag.eigenvectors() ;
     R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
     rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
     }

   kappa = R.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
   t = evec.block(  0, 2*nbas, 2*nbas, 2*nbas).adjoint()*evec.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
   nrm = std::sqrt( t.determinant()) ;

/*
  Check for convergence
*/
   R.block( 0, 0, 2*nbas, 2*nbas) = rho ;
   R.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
   R.block( 2*nbas, 0, 2*nbas, 2*nbas) = kappa.adjoint() ;
   R.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;

   std::cout << " rms difference in the density " << energy << std::endl ;
   if ( energy.real() < 1.0e-7 ) {
     std::cout << " converged in the density " << std::endl ;
     break ;
     }

 } while ( ++iter < 250) ;

  c = H_diag.eigenvectors() ;
  teig = H_diag.eigenvalues() ;
  eig = teig.real() ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }

  ring_shiekh_cg_time.end() ;

  return ;

} ;

template void ring_shiekh_cg( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::VectorXd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd,  trapezoid*&, double&, cd&, int&) ;

template < class matrix>
void tr_ring_shiekh( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> V, Eigen::Ref<Eigen::MatrixXcd> U, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, int& maxit) {

/*
  Complex Generalized Number projection with Time reversal
*/

  int iter = 1, iter_N = 0 ;
  int in, inmb = ngrid->ns() ;
  int tr_opt = 0, ci_indx ;
  double b_ll, b_ul ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> HCI_diag ;
  cd N = z0, olap ;
  cd nrm ;
  cd energy, pphase, nphase, shift ;
  cd s_theta, intE_g, intx_g, fnx, w_phi, x_phi, y_phi ;
  cd etr, gtr, dtr, d_phi ;
  cd cospi2, sinpi2 ;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::MatrixXcd R, evec, Revec ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd tmp, Y_tmp, I ;
  Eigen::MatrixXcd C_phi, G_phi, D_phi ;
  Eigen::MatrixXcd r_phi, k_phi, R_phi ;
  Eigen::MatrixXcd X_phi, Y_phi ;
  Eigen::MatrixXcd Ux, Vx, ta, tb, rho, kappa ;
  Eigen::MatrixXcd kbar_phi, t, t1 ;
  Eigen::MatrixXcd mu, mu_n, Heff, H ;
  Eigen::MatrixXcd H_CI( 2, 2), S_CI( 2, 2), CI_vec( 2, 2), CI_val( 2, 2), CI_cof( 4, 1) ;

  time_dbg ring_shiekh_projection_time = time_dbg("tr_ring_shiekh") ;

/*
  Local :
    R - Generalized density
    evec - eigenvectors/scratch space
*/

  R.resize( 4*nbas, 4*nbas) ;
  evec.resize( 4*nbas, 4*nbas) ;
  Revec.resize( 4*nbas, 4*nbas) ;
  p_rho.resize( 2*nbas, 2*nbas) ;
  p_kappa.resize( 2*nbas, 2*nbas) ;
  epsilonN.resize( 2*nbas, 2*nbas) ;
  gammaN.resize( 2*nbas, 2*nbas) ;
  lambdaN.resize( 2*nbas, 2*nbas) ;
  FockN.resize( 2*nbas, 2*nbas) ;
  DeltaN.resize( 2*nbas, 2*nbas) ;
  Y_tmp.resize( 2*nbas, 2*nbas) ;
  tmp.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  rho.resize( 2*nbas, 2*nbas) ;
  kappa.resize( 2*nbas, 2*nbas) ;
  C_phi.resize( 2*nbas, 2*nbas) ;
  G_phi.resize( 2*nbas, 2*nbas) ;
  R_phi.resize( 2*nbas, 2*nbas) ;
  D_phi.resize( 2*nbas, 2*nbas) ;
  r_phi.resize( 2*nbas, 2*nbas) ;
  k_phi.resize( 2*nbas, 2*nbas) ;
  X_phi.resize( 2*nbas, 2*nbas) ;
  Y_phi.resize( 2*nbas, 2*nbas) ;
  Ux.resize( 2*nbas, 2*nbas) ;
  Vx.resize( 2*nbas, 2*nbas) ;
  ta.resize( nbas, 2*nbas) ;
  tb.resize( nbas, 2*nbas) ;
  kbar_phi.resize( 2*nbas, 2*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  t1.resize( 2*nbas, 2*nbas) ;
  mu.resize( 4*nbas, 4*nbas) ;
  mu_n.resize( 4*nbas, 4*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;
  H.resize( 4*nbas, 4*nbas) ;

  I.setIdentity() ;
  mu.setZero() ;
  mu.block( 0, 0, 2*nbas, 2*nbas) = -I ;
  mu.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I ;
  mu_n.setIdentity() ;

/*
  We are doing time reversal projection.  Calculate cos(pi/2)
  and sin(pi/2) here one time.
  cospi2 = cos( zpi/z2) ;
  sinpi2 = sin( zpi/z2) ;
*/

  cospi2 = z0 ;
  sinpi2 = z1 ;

/*
  step 1 : Generate an initial CI vector by iterating over the operators
  we need.
*/

  for( ci_indx = 0; ci_indx < 4; ci_indx++){
    if ( ci_indx == 0){
      /* <0|P^{N}|0> */
//      std::cout << " <0|P^{N}|0> " << std::endl ;

      rho = V.conjugate()*V.transpose() ;
      kappa = V.conjugate()*U.transpose() ;
    } else if ( ci_indx == 2){
      /* <0|P^{N}T|0> */
 //     std::cout << " <0|P^{N}T|0> " << std::endl ;

      ta = V.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = V.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Vx.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Vx.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;

      rho = Vx.conjugate()*V.transpose() ;
      kappa = Vx.conjugate()*U.transpose() ;
    } else if ( ci_indx == 1){
      /* <0|TP^{N}|0> */
//      std::cout << " <0|TP^{N}|0> " << std::endl ;

      ta = V.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = V.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Vx.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Vx.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;
 
      ta = U.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = U.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Ux.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Ux.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;

      rho = V.conjugate()*Vx.transpose() ;
      kappa = V.conjugate()*Ux.transpose() ;
    } else if ( ci_indx == 3){
      /* <0|TP^{N}T|0> */
//      std::cout << " <0|TP^{N}T|0> " << std::endl ;

      ta = V.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = V.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Vx.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Vx.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;
 
      ta = U.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = U.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Ux.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Ux.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;

      rho = Vx.conjugate()*Vx.transpose() ;
      kappa = Vx.conjugate()*Ux.transpose() ;
      }

      t = kappa.inverse() ;
      evec.block( 0, 0, 2*nbas, 2*nbas) = -rho*t.conjugate() ;
      evec.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      evec.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      evec.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = t*rho ;
      olap = pfaffian_H( evec) ;
      std::cout << " olap " << olap << std::endl ;
      std::cout << " norm " << z1/olap << std::endl ;

      continue ;

/*
  Calculate the overlap and energy between the two states
*/

    intx_g = z0 ;
    intE_g = z0 ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      w_phi = ngrid->w() ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
      pphase = std::exp( z2*zi*fnx) ;
      nphase = std::exp( -z2*zi*fnx) ;

/*
  Build our rotated quantities
*/
      tmp = I + rho*( pphase - z1) ;
      C_phi = pphase*tmp.inverse() ;
      r_phi = C_phi*rho ;
      k_phi = C_phi*kappa ;
      kbar_phi = pphase*kappa*C_phi.conjugate() ;
      R_phi = I*std::exp( zi*fnx) ;
      t = R_phi*kappa ;
      t1 = t.inverse() ;
      evec.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*t1.conjugate() ;
      evec.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      evec.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
      evec.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa.inverse()*rho ;
      olap = nrm*pfaffian_H( evec) ;
      x_phi = d_phi*olap ;

      W->contract( r_phi, k_phi) ;
      evec = W->getG() ;
      G_phi = evec.block( 0, 0, 2*nbas, 2*nbas) ;
      D_phi = evec.block( 0, 2*nbas, 2*nbas, 2*nbas) ;

      t = h*r_phi ;
      etr = t.trace() ;
      t = G_phi*r_phi ;
      gtr = t.trace() ;
      t = D_phi*kbar_phi.adjoint() ;
      dtr = -t.trace()/z2 ;
      t = D_phi.transpose()*kbar_phi.conjugate() ;
      dtr -= t.trace()/z2 ;

/*
  Accumulate the projected energy and the overlap
*/

      intx_g += w_phi*x_phi ;
      intE_g += w_phi*x_phi*( etr + gtr/z2 - dtr/z2) ;
      }

    *(H_CI.data() + ci_indx) = intE_g ;
    *(S_CI.data() + ci_indx) = intx_g ;

    std::cout << " Energy " << intE_g/intx_g << std::endl ;

    ngrid->set_s() ;
 
    }

  return ;

  print_mat( H_CI, " H_CI matrix ") ;
  print_mat( S_CI, " S_CI matrix ") ;

  HCI_diag.compute( H_CI, S_CI) ;
  if ( HCI_diag.info() != Eigen::Success){
    std::cout << " Not Successful " << std::endl ;
    }

  CI_vec = HCI_diag.eigenvectors() ;
  CI_val = HCI_diag.eigenvalues() ;
  print_mat( CI_vec, " TRCI eigvec ") ;
  print_mat( CI_val, " TRCI eigval ") ;
  std::cout << " Overlap " << std::endl ;
  std::cout << CI_vec.adjoint().block( 0, 0, 1, 2)*S_CI*CI_vec.block( 0, 0, 2, 1) << std::endl ;
  std::cout << " Energy " << std::endl ;
  std::cout << CI_vec.adjoint().block( 0, 0, 1, 2)*H_CI*CI_vec.block( 0, 0, 2, 1) << std::endl ;
//  intx_g = CI_vec.adjoint().block( 0, 0, 1, 2)*S_CI*CI_vec.block( 0, 0, 2, 1) ;
//  std::cout << " Overlap " << intx_g << std::endl ;
//  intE_g = CI_vec.adjoint().block( 0, 0, 1, 2)*H_CI*CI_vec.block( 0, 0, 2, 1) ;
//  std::cout << " Energy " << intE_g << std::endl ;
//  std::cout << " <0|H|0>/<0|0> " << intE_g/intx_g << std::endl ;

  return ;

/*
  Build the prefactors for the linear combination of Time Reversal elements
*/

  CI_cof( 0, 0) = std::conj(CI_vec( 0, 0))*CI_vec( 0, 0) ;
  CI_cof( 1, 0) = std::conj(CI_vec( 1, 0))*CI_vec( 0, 0) ;
  CI_cof( 2, 0) = std::conj(CI_vec( 0, 0))*CI_vec( 1, 0) ;
  CI_cof( 3, 0) = std::conj(CI_vec( 1, 0))*CI_vec( 1, 0) ;

/*
  Loop of the Time Reversal elements and gather the heffective Hamiltonians together into one 
  MEGA EFFECTIVE HAMILTONIAN
*/

  Heff.setZero() ;


    do {

  for( ci_indx = 0; ci_indx < 4; ci_indx++){
    if ( ci_indx == 0){
      /* <0|P^{N}|0> */

      rho = V.conjugate()*V.transpose() ;
      kappa = V.conjugate()*U.transpose() ;
    } else if ( ci_indx == 2){
      /* <0|P^{N}T|0> */

      ta = V.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = V.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Vx.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Vx.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;

      rho = Vx.conjugate()*V.transpose() ;
      kappa = Vx.conjugate()*U.transpose() ;
    } else if ( ci_indx == 1){
      /* <0|TP^{N}|0> */

      ta = V.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = V.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Vx.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Vx.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;
 
      ta = U.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = U.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Ux.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Ux.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;

      rho = V.conjugate()*Vx.transpose() ;
      kappa = V.conjugate()*Ux.transpose() ;
    } else if ( ci_indx == 3){
      /* <0|TP^{N}T|0> */

      ta = V.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = V.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Vx.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Vx.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;
 
      ta = U.block( 0, 0, nbas, 2*nbas).conjugate() ;
      tb = U.block( nbas, 0, nbas, 2*nbas).conjugate() ;
      Ux.block( 0, 0, nbas, 2*nbas) = cospi2*ta + sinpi2*tb ;
      Ux.block( nbas, 0, nbas, 2*nbas) = cospi2*tb - sinpi2*ta ;

      rho = Vx.conjugate()*Vx.transpose() ;
      kappa = Vx.conjugate()*Ux.transpose() ;
      }

/*
  Loop over the number projection grid and accumulate quantities 
  Clear our accumulation matrices
*/

      p_rho = rho ;
      p_kappa = kappa ;
   
      FockN.setZero() ;
      DeltaN.setZero() ;
      Y_tmp.setZero() ;
  
      intx_g = z0 ;
      intE_g = z0 ;
  
      for ( in = 0; in < inmb; in++){
        fnx = static_cast<cd>( ngrid->q()) ;
        w_phi = ngrid->w() ;
        d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
        pphase = std::exp( z2*zi*fnx) ;
  
        R_phi = I*std::exp( zi*fnx) ;
        tmp = I + rho*( pphase - z1) ;
        C_phi = pphase*tmp.inverse() ;
  
        t = R_phi*kappa ;
        t1 = t.inverse() ;
        evec.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*t1.conjugate() ;
        evec.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
        evec.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
        evec.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa.inverse()*rho ;
        olap = nrm*pfaffian_H( evec) ;
  
        x_phi = d_phi*olap ;
        Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;
  
        intx_g += w_phi*x_phi ;
  
        }
  
      ngrid->set_s() ;
  
      Y_tmp /= intx_g ;
  
      for ( in = 0; in < inmb; in++){
        fnx = static_cast<cd>( ngrid->q()) ;
        w_phi = ngrid->w() ;
        d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
        pphase = std::exp( z2*zi*fnx) ;
        nphase = std::exp( -z2*zi*fnx) ;

/*
  Build our rotated quantities
*/

        tmp = I + rho*( pphase - z1) ;
        C_phi = pphase*tmp.inverse() ;
        r_phi = C_phi*rho ;
        k_phi = C_phi*kappa ;
        kbar_phi = pphase*kappa*C_phi.conjugate() ;
        R_phi = I*std::exp( zi*fnx) ;
        t = R_phi*kappa ;
        t1 = t.inverse() ;
        evec.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*t1.conjugate() ;
        evec.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
        evec.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
        evec.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa.inverse()*rho ;
        olap = nrm*pfaffian_H( evec) ;
        x_phi = d_phi*olap ;
        y_phi = x_phi/intx_g ;
        Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi - Y_tmp ;
  
        W->contract( r_phi, k_phi) ;
        evec = W->getG() ;
        G_phi = evec.block( 0, 0, 2*nbas, 2*nbas) ;
        D_phi = evec.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
  
        t = h*r_phi ;
        etr = t.trace() ;
        epsilonN = w_phi*y_phi*(Y_phi*etr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi) ;
        t = G_phi*r_phi ;
        gtr = t.trace() ;
        gammaN = w_phi*y_phi*(Y_phi*gtr/z2 + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*G_phi*C_phi) ;
        t = D_phi*kbar_phi.adjoint() ;
        dtr = -t.trace()/z2 ;
        t = D_phi.transpose()*kbar_phi.conjugate() ;
        dtr -= t.trace()/z2 ;
        lambdaN = -w_phi*y_phi*(Y_phi*dtr/z2 + z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint()) ;
  
        FockN += epsilonN + gammaN + lambdaN ;
        DeltaN += w_phi*nphase*y_phi*C_phi*D_phi ;

/*
  Accumulate the overlap and projected energy
*/
        }

      ngrid->set_s() ;
  
      t = (FockN + FockN.adjoint())/z2 ;
      FockN = t ;
      t = (DeltaN - DeltaN.transpose())/z2 ;
      DeltaN = t ;

      Heff.block( 0, 0, 2*nbas, 2*nbas) += CI_cof( ci_indx, 0)*FockN ;
      Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) += CI_cof( ci_indx, 0)*DeltaN ;
      Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) += -CI_cof( ci_indx, 0)*FockN.conjugate() ;
      Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) += CI_cof( ci_indx, 0)*DeltaN.adjoint() ;

    }

    if ( iter < 5 ){
/*
  Slowly scale down the level shifting
*/
      shift = lshift/static_cast<cd>(iter) ;
      }

    Heff += -shift*R ;

/*
  Adjust the chemical potential until the density gives the proper number
  of particles.
*/

    if ( true ) {
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
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
        rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
        N = rho.trace() ;
        } while ( static_cast<double>(N.real()) > nele) ;
    
      std::cout << " Lower limit chemical potential " << b_ll << std::endl ;
      iter_N = 0 ;
    
      do {
        b_ul += d3 ;
        H = Heff + static_cast<cd>(b_ul)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
        rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
        N = rho.trace() ;
        } while ( static_cast<double>(N.real()) < nele) ;
    
      std::cout << " upper limit chemical potential " << b_ul << std::endl ;
      iter_N = 0 ;
    
      lambda = ( b_ll + b_ul)/d2 ;
    
      while ( iter_N++ < 30) {
        H = Heff + static_cast<cd>(lambda)*mu ;
        H_diag.compute( H) ;
        evec = H_diag.eigenvectors() ;
        R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
        rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
        N = rho.trace() ;
        if ( std::abs(static_cast<double>(N.real()) - nele) < 1.0e-5){
          std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
          std::cout << "    chemical potential: " << lambda << std::endl ;
          std::cout << "    Particle Number : " << N.real() << std::endl << std::endl ;
          break ;
        } else {
  
          if ( static_cast<double>(N.real()) - nele < d0){
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
        }
   } else {
/*
  If there is no chemical potential then simply diagonalize the Hamiltonian
*/
     H_diag.compute( Heff) ;
     evec = H_diag.eigenvectors() ;
     R = evec.block( 0, 0, 4*nbas, 2*nbas)*evec.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;
     rho = R.block( 0, 0, 2*nbas, 2*nbas) ;
     }

   kappa = R.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
   t = evec.block(  0, 2*nbas, 2*nbas, 2*nbas).adjoint()*evec.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
   nrm = std::sqrt( t.determinant()) ;

/*
  Check for convergence
*/

   R.block( 0, 0, 2*nbas, 2*nbas) = rho ;
   R.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
   R.block( 2*nbas, 0, 2*nbas, 2*nbas) = kappa.adjoint() ;
   R.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;

   std::cout << " rms difference in the density " << energy << std::endl ;
   if ( energy.real() < 1.0e-7 ) {
     std::cout << " converged in the density " << std::endl ;
     break ;
     }

/*
  Get the new U and V and continue the cycle.
*/
  V = evec.block( 0, 0, 2*nbas, 2*nbas) ;
  U = evec.block( 2*nbas, 0, 2*nbas, 2*nbas) ;

 } while ( ++iter < 250) ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }

  ring_shiekh_projection_time.end() ;

  return ;

} ;

template void tr_ring_shiekh( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd,  trapezoid*&, double&, int&) ;

//void rrHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, trapezoid*& ngrid, int& maxit, Eigen::Ref<Eigen::MatrixXcd> xs, Eigen::Ref<Eigen::MatrixXcd> xsi) {
//
///*
//  Various debugging for the number projected HFB.
//*/
//
//  int iter = 0 ; //, iter_N = 0 ;
//  int in ; //, diis_on = 0 ;
//  int inmb = ngrid->ns() ;
////  diis<cd> W = diis<cd>( 16*nbas*nbas, 5) ;
////  double b_ll, b_ul, nele = d4 ;
//  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
//  Eigen::BDCSVD<Eigen::MatrixXcd> H_SVD ;
//  cd lshift = -z4, x ;
////  cd N = z0 ;
//  cd z2pi ( d2*pi, d0) ;
//  cd vac_nrm, energy = z0 ; /* This is the vacuum normalization */
//  cd s_theta, intE_g, intx_g, x_g, p_energy ;
//  cd blah3, blah2, blah, fnx, e_theta, fnw;
//  std::vector<cd> cnv_den, cnv_energy ;
//  Eigen::MatrixXcd t ;
//  Eigen::MatrixXcd nX ;
//  Eigen::MatrixXcd Rp ;
//  Eigen::MatrixXcd Heff, H ;
//  Eigen::MatrixXcd Sdotrho, Sdotkappa ;
//  Eigen::MatrixXcd Hdotrho, Hdotkappa ;
//  Eigen::MatrixXcd f_theta ;
//  Eigen::MatrixXcd kappa_bar ;
//  Eigen::MatrixXcd p_rho, p_kappa ;
//  Eigen::MatrixXcd Rrho ;
//  Eigen::MatrixXcd Rkappa ;
//  Eigen::MatrixXcd rho_i, kappa_i ;
//  Eigen::MatrixXcd Rrho_i ;
//  Eigen::MatrixXcd Rkappa_i ;
//  Eigen::MatrixXcd tmp ;
//  Eigen::MatrixXcd r_theta ;
//  Eigen::MatrixXcd k_theta ;
//  Eigen::MatrixXcd mu, mu_n ;
//  Eigen::MatrixXcd I, G, D, D_bar, Vv, lvlshft ;
//  Eigen::MatrixXcd Zm1, M11, M22 ;
//  Eigen::IOFormat tmpfmt(4, 0, ", ", "\n", "[", "]") ;
//
///*
//  Eigen Vectors to accumulate quantities
//*/
//
//  Eigen::VectorXcd OHg( inmb), OSg( inmb), d_g( inmb) ;
//  Eigen::MatrixXcd C_theta ;
//  Eigen::MatrixXcd C_theta_i ;
//  time_dbg cgHFB_projection_time = time_dbg("rrHFB_projection") ;
//  t.resize( 2*nbas, 2*nbas) ;
//  h.resize( 2*nbas, 2*nbas) ;
//  f_theta.resize( 2*nbas, 2*nbas) ;
//  r_theta.resize( 2*nbas, 2*nbas) ;
//  k_theta.resize( 2*nbas, 2*nbas) ;
//  Zm1.resize( 2*nbas, 2*nbas) ;
//  mu.resize( 4*nbas, 4*nbas) ;
//  mu_n.resize( 4*nbas, 4*nbas) ;
//  Vv.resize( 4*nbas, 4*nbas) ;
//  tmp.resize( 2*nbas, 2*nbas) ;
///* Normal rho and kappa */
//  p_rho.resize( 2*nbas, 2*nbas) ;
//  p_kappa.resize( 2*nbas, 2*nbas) ;
///* Rotated rho and kappa */
//  Rkappa.resize( 2*nbas, 2*nbas) ;
//  Rrho.resize( 2*nbas, 2*nbas) ;
///* Inverse Rho and kappa  */
//  kappa_i.resize( 2*nbas, 2*nbas) ;
//  rho_i.resize( 2*nbas, 2*nbas) ;
///* Rotated Inverse  */
//  Rkappa_i.resize( 2*nbas, 2*nbas) ;
//  Rrho_i.resize( 2*nbas, 2*nbas) ;
////
//  C_theta.resize( 2*nbas, 2*nbas) ;
//  C_theta_i.resize( 2*nbas, 2*nbas) ;
//  Sdotrho.resize( 2*nbas, 2*nbas) ;
//  Sdotkappa.resize( 2*nbas, 2*nbas) ;
//  Hdotrho.resize( 2*nbas, 2*nbas) ;
//  Hdotkappa.resize( 2*nbas, 2*nbas) ;
//  Rp.resize( 4*nbas, 4*nbas) ;
//  Heff.resize( 4*nbas, 4*nbas) ;
//  H.resize( 4*nbas, 4*nbas) ;
//  nX.resize( 2*nbas, 2*nbas) ;
//  I.resize( 2*nbas, 2*nbas) ;
//  I.setIdentity() ;
//  Heff.setZero() ;
//  Rp.setZero() ;
//  G.resize( 2*nbas, 2*nbas) ;
//  D.resize( 2*nbas, 2*nbas) ;
//  D_bar.resize( 2*nbas, 2*nbas) ;
//  M11.resize( 2*nbas, 2*nbas) ;
//  M22.resize( 2*nbas, 2*nbas) ;
//  lvlshft.resize( 4*nbas, 4*nbas) ;
//  mu.setZero() ;
//  mu.block( 0, 0, 2*nbas, 2*nbas) = -I ;
//  mu.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I ;
//  mu_n.setIdentity() ;
//
//  vac_nrm = z1 ;
//  cd sgn = std::pow( -z1, 4*nbas*( 4*nbas + 1)/2) ;
//
//  print_mat( rho, " initial rho") ;
//  print_mat( kappa, " initial kappa") ;
//
//  /*
//    Put into the orthonormal basis and prepare level shifting density
//  */
//
//  transform( 2, xsi, rho) ;
//  transform( 2, xsi, kappa) ;
//  lvlshft.block( 0, 0, 2*nbas, 2*nbas) = rho ;
//  lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
//  lvlshft.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
//  lvlshft.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
//  p_rho = rho ;
//  p_kappa = kappa ;
//
//  do {
//    /*
//      Loop over the number projection grid and accumulate quantities 
//    */
//    rho = p_rho ;
//    rho_i = rho.inverse() ;
//    kappa_i = kappa.inverse() ;
//
//    /*
//      Clear our accumulation matrices
//    */
//
//    Sdotrho.setZero() ;
//    Sdotkappa.setZero() ;
//    Hdotrho.setZero() ;
//    Hdotkappa.setZero() ;
//
//    intx_g = z0 ;
//    intE_g = z0 ;
//
//    for ( in = 0; in < inmb; in++){
//      fnx = static_cast<cd>( ngrid->q()) ;
//      d_g(in) = ngrid->w()*std::exp( -zi*fnx*z4) ;
//
//      /*
//        Build our rotated quantities
//      */
//
//      nX = I*std::exp( zi*fnx) ;
//      tmp = nX.inverse() ;
//      Rrho = nX*rho ;
//      Rrho_i = Rrho.inverse() ;
//      Rkappa = nX*kappa ;
//      Rkappa_i = Rkappa.inverse() ;
//      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
//      C_theta = C_theta_i.inverse() ;
//      r_theta = Rrho*C_theta*rho ;
//      k_theta = Rrho*C_theta*kappa ;
//      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;
//
//      /*
//        Get the energy expression.
//        Start by putting into the non orthogonal basis
//      */
//
//      transform( 2, xs, r_theta) ;
//      transform( 2, xs, k_theta) ;
//      transform( 2, xs, kappa_bar) ;
//
//      ctr2eg( intarr, r_theta, G, nbas) ;
//      ctrPairg( intarr, k_theta, D, nbas) ;
//      ctrPairg( intarr, kappa_bar, D_bar, nbas) ;
//      D /= z2 ;
//      D_bar /= z2 ;
//      transform( 2, xsi, r_theta) ;
//      transform( 2, xsi, k_theta) ;
//      transform( 2, xsi, kappa_bar) ;
//      transform( 2, xs, G) ;
//      transform( 2, xs, D) ;
//      transform( 2, xs, D_bar) ;
//      f_theta = h + G ;
//      t = (h + f_theta)*r_theta ;
//      OHg( in) = t.trace()/z2 ;
//      t = kappa_bar*D/z2 ;
//      OHg( in) += t.trace() ;
//
//      /* 
//        Get the overlap and the necessary matrix inverses 
//      */
//
//      Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
//      Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
//      Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
//      Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
////      tmp = -Rrho*Rkappa_i.conjugate() + rho_i*kappa_i ;
////      M11 = tmp.inverse() ;
////      M22 = rho_i*kappa*(I - M11*rho_i*kappa) ;
////      Zm1 = -Rkappa_i.conjugate()*M11*nX + M22*kappa_i ;
//      Zm1 = nX*rho*C_theta + C_theta*rho*nX ;
// 
//      if ( in == 0 ) {
//        vac_nrm = std::abs(z1/pfaffian_H( Rp)) ;
//        std::cout << " vac_nrm at in = 0 " << std::endl ;
//        std::cout << vac_nrm << std::endl ;
//        Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
//        Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
//        Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
//        Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
//        }
//
//      OSg( in) = vac_nrm*sgn*pfaffian_H( Rp) ;
//
//      /*
//        Sdotrho
//      */
//      Sdotrho += d_g( in)*OSg( in)*Zm1/z2 ;
//
//      /*
//        Hdotrho
//      */
//      tmp = C_theta*rho*f_theta*nX ;
//      tmp += f_theta*Rrho*C_theta ;
//      tmp += -r_theta*f_theta*Rrho*C_theta ;
//      tmp += -C_theta*rho*f_theta*r_theta*nX ;
//
//      tmp += r_theta*D*Rkappa.conjugate()*C_theta/z2 ;
//      tmp += -D*Rkappa.conjugate()*C_theta/z2 ;
//      tmp += -C_theta*rho*D*kappa_bar*nX/z2 ;
//
//      tmp += -C_theta*kappa*D_bar*nX/z2 ;
//      tmp += k_theta*D_bar*Rrho*C_theta/z2 ;
//      tmp += C_theta*kappa*D_bar*r_theta*nX/z2 ;
//
//      Hdotrho += d_g( in)*OSg( in)*(OHg( in)*Zm1/z2 + tmp) ;
///*  
//  Sdotkappa
//*/
//      Zm1 = C_theta*kappa*nX.conjugate() ;
////      print_mat( Zm1, " kappa det derivative") ;
////      Zm1 = Rkappa_i.conjugate()*M11*Rrho*kappa_i.conjugate() ;
////      print_mat( Zm1, " kappa derivative S") ;
//      Sdotkappa += d_g( in)*OSg( in)*Zm1/z2 ;
///*  
//  Hdotkappa
//*/
//      tmp = C_theta*rho*f_theta*k_theta*nX.conjugate() ;
//      tmp += C_theta*kappa*D_bar*k_theta*nX.conjugate()/z2 ;
//      tmp += -C_theta*rho*D*nX.conjugate()/z2 ;
//      tmp += -C_theta*rho*D*Rkappa.conjugate()*C_theta*kappa*nX.conjugate()/z2 ;
//
//      Hdotkappa += d_g( in)*OSg( in)*(OHg( in)*Zm1/z2 + tmp) ;
///*   
//  Accumulate the overlap and projected energy
//*/  
//      intx_g += OSg( in)*d_g( in) ;
//      intE_g += OSg( in)*OHg( in)*d_g( in) ;
//      }
//
//    ngrid->set_s() ;
//    energy = intE_g/intx_g ;
//    cnv_energy.push_back( energy) ;
//    tmp = Hdotrho ;
//    Hdotrho = (tmp + tmp.adjoint())/z2 ;
//    tmp = Sdotrho ;
//    Sdotrho = (tmp + tmp.adjoint())/z2 ;
////    tmp = Hdotrho - intE_g*Sdotrho/intx_g ;
////    tmp /= intx_g ;
////    print_mat( tmp, "Hdotrho noao") ;
////    transform( 2, xsi, tmp) ;
////    print_mat( tmp, "Hdotrho oao") ;
//    tmp = Hdotkappa ;
//    Hdotkappa = (tmp - tmp.transpose())/z2 ;
//    tmp = Sdotkappa ;
//    Sdotkappa = (tmp - tmp.transpose())/z2 ;
////    tmp = Hdotkappa - intE_g*Sdotkappa/intx_g ;
////    tmp /= intx_g ;
////    print_mat( tmp, "Hdotkappa noao") ;
////    transform( 2, xsi, tmp) ;
////    print_mat( tmp, "Hdotkappa oao") ;
//
//    Heff.block( 0, 0, 2*nbas, 2*nbas) = Hdotrho/intx_g - intE_g*Sdotrho/std::pow( intx_g, z2) ;
//    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = Hdotkappa/intx_g - intE_g*Sdotkappa/std::pow( intx_g, z2) ;
//    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
//    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
//    print_mat( Heff, " Heff ") ;
//    mu_n.setIdentity() ;
//
//    Heff += lshift*lvlshft ;
//    Heff -= lshift*(mu_n - lvlshft) ;
//    /*
//      Adjust the chemical potential until the density gives the proper number
//      of particles.
//    */
//
//    iter_N = 0 ;
//
//    b_ul = lambda ;
//    b_ll = lambda ;
//
//    /*
//       Set some initial limits
//    */
//
//    do {
//      b_ll -= d3 ;
//      H = Heff + static_cast<cd>(b_ll)*mu ;
//      H_diag.compute( H) ;
//      Vv = H_diag.eigenvectors() ;
//      rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//      N = rho.trace() ;
//      } while ( static_cast<double>(N.real()) > static_cast<double>(nele)) ;
//  
//    std::cout << " Lower limit chemical potential " << b_ll << std::endl ;
//  
//    do {
//      b_ul += d3 ;
//      H = Heff + static_cast<cd>(b_ul)*mu ;
//      H_diag.compute( H) ;
//      Vv = H_diag.eigenvectors() ;
//      rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//      N = rho.trace() ;
//      } while ( static_cast<double>(N.real()) < static_cast<double>(nele)) ;
//  
//    std::cout << " upper limit chemical potential " << b_ul << std::endl ;
//  
//    lambda = ( b_ll + b_ul)/d2 ;
//  
//    while ( iter_N++ < 30) {
//      H = Heff + static_cast<cd>(lambda)*mu ;
//      H_diag.compute( Heff) ;
//      Vv = H_diag.eigenvectors() ;
//      rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//      N = rho.trace() ;
//      std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
//      std::cout << "    chemical potential: " << lambda << std::endl ;
//      std::cout << "    Particle Number : " << N.real() << std::endl << std::endl ;
//      if ( std::abs(static_cast<double>(N.real()) - static_cast<double>(nele)) < 1.0e-5){
//        break ;
//      } else {
//  
/*
  Bisection method
*/

//        if ( static_cast<double>(N.real()) - static_cast<double>(nele) < d0){
/*
  Too few electrons. Increase the chemical potential
*/
//
//          b_ll = lambda ;
//          lambda = (b_ul + b_ll)/d2 ;
//        } else {
//          b_ul = lambda ;
//          lambda = (b_ul + b_ll)/d2 ;
//          }
//        }
//      }

//   H_diag.compute( Heff) ;
//   std::cout << " H_diag.eigenvalues() " << std::endl ;
//   std::cout << H_diag.eigenvalues() << std::endl ;
//   Vv = H_diag.eigenvectors() ;
//   rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//
//   kappa = Vv.block(  2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//   tmp = rho ;
//   rho = (tmp + tmp.adjoint())/z2 ;
//   tmp = kappa ;
//   kappa = (tmp - tmp.transpose())/z2 ;
//   lvlshft.block( 0, 0, 2*nbas, 2*nbas) = rho ;
//   lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
//   lvlshft.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
//   lvlshft.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
//   if( iter > 5 && diis_on == 0){
//     std::cout << "diis is on" << std::endl ;
//     diis_on = 1 ;
//     }
//   W.update( lvlshft, diis_on) ;
//   rho = lvlshft.block( 0, 0, 2*nbas, 2*nbas) ;
//   kappa = lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) ;

/*
  Check for convergence
*/
//   t = rho - p_rho ;
//   energy = t.norm() ;
//   t = kappa - p_kappa ;
//   energy += t.norm() ;
//   cnv_den.push_back( energy) ;
//
//   std::cout << " rms difference in the density " << energy << std::endl ;
//   if ( energy.real() < 1.0e-7 ) {
//     std::cout << " converged in the density " << std::endl ;
//     break ;
//     }
//
//   p_rho = rho ;
//   p_kappa = kappa ;
//
//   print_mat( rho, " rho ") ;
//   print_mat( kappa, " kappa ") ;
// 
// } while ( ++iter < 1) ;
//
//  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
//    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
//    }
//
//  cgHFB_projection_time.end() ;
//
//  return ;
//
//} ;

//void cgHFB_projection_dbg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, trapezoid*& ngrid, int& maxit) {
//
///*
//  Various debugging for the number projected HFB.
//*/
//
//  int iter = 0 ;
//  int in ;
//  int inmb = ngrid->ns() ;
//  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
//  cd z2pi ( d2*pi, d0) ;
//  cd vac_nrm, energy = z0 ; /* This is the vacuum normalization */
//  cd blah3, blah2 = z0, blah = z0, e_theta, fnw, fnx, s_theta, intE_g = z0, intx_g = z0, x_g, p_energy=z0 ;
//  std::vector<cd> cnv_den ;
//  Eigen::MatrixXcd t ;
//  Eigen::MatrixXcd nX ;
//  Eigen::MatrixXcd Rp ;
//  Eigen::MatrixXcd Heff ;
//  Eigen::MatrixXcd Sdotrho ;
//  Eigen::MatrixXcd Sdotkappa ;
//  Eigen::MatrixXcd Hdotrho ;
//  Eigen::MatrixXcd Hdotkappa ;
//  Eigen::MatrixXcd f_theta ;
//  Eigen::MatrixXcd kappa_bar ;
//  Eigen::MatrixXcd p_rho ;
//  Eigen::MatrixXcd p_kappa ;
//  Eigen::MatrixXcd Rrho ;
//  Eigen::MatrixXcd Rkappa ;
//  Eigen::MatrixXcd rho_i ;
//  Eigen::MatrixXcd kappa_i ;
//  Eigen::MatrixXcd Rrho_i ;
//  Eigen::MatrixXcd Rkappa_i ;
//  Eigen::MatrixXcd tmp ;
//  Eigen::MatrixXcd Rp_i ;
//  Eigen::MatrixXcd r_theta ;
//  Eigen::MatrixXcd k_theta ;
//  Eigen::MatrixXcd I, G, D, D_bar, Vv ;
//  Eigen::MatrixXcd Zm1, M11, M22, lvlshft ;
//
///*
//  Eigen Vectors to accumulate quantities
//*/
//
//  Eigen::VectorXcd OHg( inmb), OSg( inmb), d_g( inmb) ;
//
///*
//  C_theta = [-kappa*R^{*}kappa^{*} + rho*R*rho]^{-1}
//*/
//
//  Eigen::MatrixXcd C_theta ;
//  Eigen::MatrixXcd C_theta_i ;
//  time_dbg cgHFB_projection_time = time_dbg("cgHFB_projection") ;
//
//  Rp_i.resize( 4*nbas, 4*nbas) ;
//  t.resize( 2*nbas, 2*nbas) ;
//  f_theta.resize( 2*nbas, 2*nbas) ;
//  r_theta.resize( 2*nbas, 2*nbas) ;
//  k_theta.resize( 2*nbas, 2*nbas) ;
//  Zm1.resize( 2*nbas, 2*nbas) ;
//  Vv.resize( 4*nbas, 4*nbas) ;
//  tmp.resize( 2*nbas, 2*nbas) ;
///* Normal rho and kappa */
//  p_rho.resize( 2*nbas, 2*nbas) ;
//  p_kappa.resize( 2*nbas, 2*nbas) ;
///* Rotated rho and kappa */
//  Rkappa.resize( 2*nbas, 2*nbas) ;
//  Rrho.resize( 2*nbas, 2*nbas) ;
///* Inverse Rho and kappa  */
//  kappa_i.resize( 2*nbas, 2*nbas) ;
//  rho_i.resize( 2*nbas, 2*nbas) ;
///* Rotated Inverse  */
//  Rkappa_i.resize( 2*nbas, 2*nbas) ;
//  Rrho_i.resize( 2*nbas, 2*nbas) ;
//
//  C_theta.resize( 2*nbas, 2*nbas) ;
//  C_theta_i.resize( 2*nbas, 2*nbas) ;
//  Sdotrho.resize( 2*nbas, 2*nbas) ;
//  Sdotkappa.resize( 2*nbas, 2*nbas) ;
//  Hdotrho.resize( 2*nbas, 2*nbas) ;
//  Hdotkappa.resize( 2*nbas, 2*nbas) ;
//  Rp.resize( 4*nbas, 4*nbas) ;
//  Heff.resize( 4*nbas, 4*nbas) ;
//  nX.resize( 2*nbas, 2*nbas) ;
//  I.resize( 2*nbas, 2*nbas) ;
//  I.setIdentity() ;
//  Heff.setZero() ;
//  Rp.setZero() ;
//  G.resize( 2*nbas, 2*nbas) ;
//  D.resize( 2*nbas, 2*nbas) ;
//  D_bar.resize( 2*nbas, 2*nbas) ;
//  M11.resize( 2*nbas, 2*nbas) ;
//  M22.resize( 2*nbas, 2*nbas) ;
//  lvlshft.resize( 4*nbas, 4*nbas) ;
//
///*
//  I guess we will go on RMS density for convergence.
//*/
//  
////  print_mat( rho, "rho") ;
////  print_mat( kappa, " kappa ") ;
//
///*
//  Get the proper sign for using the pfaffian.
//
//  Vacuum normalization for quasi-particles is given by the expression sqrt( det(U^{t}U).
//  Calculate it here for the initial iteration.
//
//*/
//  cd sgn = std::pow( -z1, 4*nbas*( 4*nbas + 1)/2) ;
////  std::cout << " Particle Number " << rho.trace() << std::endl ;
///*
//  tmp = w.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint()*w.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
//  vac_nrm = std::sqrt( tmp.determinant()) ;
//  std::cout << " vac_nrm " << std::endl ;
//  std::cout << vac_nrm << std::endl ;
//  since we are populating from a HF determinant, we have no initial U to get the vacuum normalization
//*/
//
//  vac_nrm = z1 ; 
//
//  do {
//    /*
//      Loop over the number projection grid and accumulate quantities
//    */
//    p_energy = energy ;
//    p_rho = rho ;
//    p_kappa = kappa ;
//    rho_i = rho.inverse() ;
//    kappa_i = kappa.inverse() ;
//    /*
//      Clear our accumulation matrices
//    */
//    Sdotrho.setZero() ;
//    Sdotkappa.setZero() ;
//    Hdotrho.setZero() ;
//    Hdotkappa.setZero() ;
//    intx_g = z0 ;
//    intE_g = z0 ;
//    for ( in = 0; in < inmb; in++){
//      fnx = static_cast<cd>( ngrid->q()) ;
//      d_g(in) = ngrid->w()*std::exp( -zi*fnx*z4)/z2pi ;
//
///*
//  Build our rotated quantities
//*/
//      nX = I*std::exp( zi*fnx) ;
//      Rrho = nX*rho ;
//      Rrho_i = Rrho.inverse() ;
//      Rkappa = nX*kappa ;
//      Rkappa_i = Rkappa.inverse() ;
//      C_theta_i = rho*Rrho - kappa*Rkappa.conjugate() ;
//      C_theta = C_theta_i.inverse() ;
//      r_theta = Rrho*C_theta*rho ;
//      k_theta = Rrho*C_theta*kappa ;
//      kappa_bar = -Rkappa.conjugate()*C_theta*rho ;
////
////      print_mat( nX, " Rotation matrix ") ;
////      print_mat( Rrho, " Rrho ") ;
////      print_mat( Rkappa, " Rkappa ") ;
////      print_mat( Rrho_i, " Rrho_i ") ;
////      print_mat( Rkappa_i, " Rkappa_i ") ;
////      print_mat( C_theta_i, " C_theta_i ") ;
////      std::cout << " det(C_theta_i) " << std::endl ;
////      std::cout << C_theta_i.determinant() << std::endl ;
////      print_mat( C_theta, " C_theta ") ; 
////      std::cout << " r_theta " << std::endl ;
////      print_mat( r_theta) ;
////      std::cout << " k_theta " << std::endl ;
////      print_mat( k_theta) ;
////      std::cout << " kappa_bar " << std::endl ;
////      print_mat( kappa_bar) ;
////
///* 
//  Get the energy expression
//*/
//    ctr2eg( intarr, r_theta, G, nbas) ;
////    print_mat( G, " G") ;
//    ctrPairg( intarr, k_theta, D, nbas) ;
////    print_mat( D, " D") ;
//    ctrPairg( intarr, kappa_bar, D_bar, nbas) ;
////    print_mat( D_bar, " D_bar") ;
//    f_theta = h + G/z2 ;
////    print_mat( f_theta. " f_theta") ;
//    t = f_theta*r_theta ;
//    OHg( in) = t.trace() ;
//    t = D*kappa_bar ;
//    OHg( in) += t.trace()/d4 ;
////    std::cout << " H " << std::endl ;
////    std::cout << OHg( in) << std::endl ;
//    f_theta = h + G ;
///* 
//  Get the overlap and the necessary matrix inverses 
//*/
//    Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
//    Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
//    Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
//    Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
//    tmp = -Rrho*Rkappa_i.conjugate() + rho_i*kappa ;
//    M11 = tmp.inverse() ;
//    M22 = rho_i*kappa*(I - M11*rho_i*kappa) ;
//    Zm1 = -Rkappa_i.conjugate()*M11*nX + M22*kappa_i ;
//    OSg( in) = vac_nrm*sgn*pfaffian_H( Rp) ;
///*
//  Sdotrho
//*/
//    Sdotrho += d_g( in)*OSg( in)*Zm1/z2 ;
///*
//  Hdotrho
//*/
//    tmp = C_theta*rho*f_theta*nX ;
//    tmp += f_theta*Rrho*C_theta ;
//    tmp += -r_theta*f_theta*Rrho*C_theta ;
//    tmp += -C_theta*rho*f_theta*r_theta*nX ;
//    tmp += r_theta*D*Rkappa.conjugate()*C_theta/z4 ;
//    tmp += -D*Rkappa.conjugate()*C_theta/z4 ;
//    tmp += -C_theta*rho*D*kappa_bar*nX/z4 ;
////
//    tmp += C_theta*kappa*D_bar*nX/z4 ;
//    tmp += -k_theta*D_bar*Rrho*C_theta/z4 ;
//    tmp += -C_theta*kappa*D_bar*r_theta*nX/z4 ;
//    Hdotrho += d_g( in)*OSg( in)*(OHg( in)*Zm1/z2 + tmp) ;
///*
//  Sdotkappa
//*/
//
//    Zm1 = Rkappa_i.conjugate()*M11*Rrho*kappa_i.conjugate() ;
//    Sdotkappa += d_g( in)*OSg( in)*(Zm1 - Zm1.transpose())/z4 ;
///*
//  Hdotkappa
//*/
//    tmp = C_theta*rho*f_theta*k_theta*nX.conjugate() ;
//    tmp += C_theta*kappa*D_bar*k_theta*nX.conjugate()/z4 ;
//    tmp += -C_theta*rho*D*nX.conjugate()/z4 ;
//    tmp += -C_theta*rho*D*Rkappa.conjugate()*C_theta*kappa*nX.conjugate()/z4 ;
//    Hdotkappa += d_g( in)*OSg( in)*(OHg( in)*(Zm1 - Zm1.transpose())/z4 + (tmp - tmp.transpose())/z2) ;
///* 
//  Accumulate the overlap and projected energy
//*/
//    intx_g += OSg( in)*d_g( in) ;
//    intE_g += OSg( in)*OHg( in)*d_g( in) ;
//    }
//
//    ngrid->set_s() ;
//
//    tmp = Hdotrho/intx_g  - intE_g*Sdotrho/std::pow( intx_g, z2) ;
//    Heff.block( 0, 0, 2*nbas, 2*nbas) = (tmp + tmp.adjoint())/z2 ;
//    tmp = Hdotkappa/intx_g - intE_g*Sdotkappa/std::pow( intx_g, z2) ;
//    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = (tmp - tmp.transpose())/z2 ;
//    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;
//    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
////    print_mat( Heff, " Heff") ;
///*
//  Level shifting because how else will this work
//*/
//
//    lvlshft.block( 0, 0, 2*nbas, 2*nbas) = rho ;
//    lvlshft.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
//    lvlshft.block( 2*nbas, 0, 2*nbas, 2*nbas) = -kappa.conjugate() ;
//    lvlshft.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;
//    Heff += lvlshft/z4 ;
//
//    H_diag.compute( Heff) ;
//    Vv = H_diag.eigenvectors() ;
//    print_mat( Vv, " eigenvectors") ;
//    std::cout << "Are the vectors orthonormal" << std::endl ;
//    std::cout << Vv.adjoint()*Vv << std::endl ;
//    rho = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//    tmp = (rho + rho.adjoint())/z2 ;
//    rho = tmp ;
////    print_mat( rho, " rho") ;
//    kappa = Vv.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas).conjugate()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose() ;
//    tmp = (kappa - kappa.transpose())/z2 ;
//    kappa = tmp ;
//    tmp = Vv.block( 0, 2*nbas, 2*nbas, 2*nbas).adjoint()*Vv.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
//    vac_nrm = std::sqrt( tmp.determinant()) ;
////    std::cout << " vac_nrm " << std::endl ;
////    std::cout << vac_nrm << std::endl ;
////    print_mat( kappa, " kappa projected") ;
///*
//  Check for convergence
//*/  
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
//
//    ngrid->set_s() ;
//
///*
//  Energy
//*/
//   energy = energy/intx_g ;
//   std::cout << "PHF Energy" << std::endl ;
//   std::cout << energy << std::endl ;
//
//   t = rho - p_rho ;
//   energy = ( t*t.adjoint()).norm() ;
//   t = kappa - p_kappa ;
//   energy += ( t*t.adjoint()).norm() ;
//   cnv_den.push_back( energy) ;
//   if ( energy.real() < 0.00001 ) { break ;}
//   } while ( ++iter < 1) ;
//
//  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
//    std::cout << cnv_den[xq] << std::endl ;
//    }
//
//  cgHFB_projection_time.end() ;
//
//  return ;
//
//} ;

//cd pf_overlap( int& nbas, cd& c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, Eigen::Ref<Eigen::MatrixXcd> R){
//
///*
//  nbas - number of basis functions
//  c - normalization term
//  rho - unrotated density
//  kappa - unrotated abnormal density
//  R - rotation matrix
//*/
//
//  cd pf ;
//  Eigen::MatrixXcd kappa_i ;
//  Eigen::MatrixXcd Rrho ;
//  Eigen::MatrixXcd Rkappa ;
//  Eigen::MatrixXcd Rkappa_i ;
//  Eigen::MatrixXcd Rp ;
//  Eigen::MatrixXcd I ;
//
//  Rp.resize( 4*nbas, 4*nbas) ;
//  kappa_i.resize( 2*nbas, 2*nbas) ;
//  Rrho.resize( 2*nbas, 2*nbas) ;
//  Rkappa.resize( 2*nbas, 2*nbas) ;
//  Rkappa_i.resize( 2*nbas, 2*nbas) ;
//  I.resize( 2*nbas, 2*nbas) ;
//  I.setIdentity() ;
//  Rrho = R*rho ;
//  kappa_i = kappa.inverse() ;
//  Rkappa = R*kappa ;
//  Rkappa_i = Rkappa.inverse() ;
//
//  Rp.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
//  Rp.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
//  Rp.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
//  Rp.block( 0, 0, 2*nbas, 2*nbas) = -Rrho*Rkappa_i.conjugate() ;
//
//  pf = pfaffian_A( Rp) ;
//
//  Rkappa_i.resize( 0, 0) ;
//  Rkappa.resize( 0, 0) ;
//  Rrho.resize( 0, 0) ;
//  kappa_i.resize( 0, 0) ;
//  Rp.resize( 0, 0) ;
//
//  return pf ;
//
//} ;
//

void jorge_guess( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, double N, double& norm) {
  int i = 0, j, nbas =  p.rows() ;
  double s1 = d1, s2, s3, n1, n2 ;
  Eigen::VectorXd u( nbas), v( nbas), v1( nbas), v2( nbas) ;
  Eigen::MatrixXd t( nbas, nbas) ;
  std::ifstream F_IN ;
  std::string line ;
/*
  This is a very specialized routine to read in a guess generated from 
  a file generated by Jorge's PBCS code.

  Read in a list of v values from a PBCS calculation.  We don't need u
  because u**2 + v**2 = 1.

  Then we do as Tom Henderson does and scale the v until v.dot(v) = # particles
  Once this is completed, the u's are generated and the normal and pairing
  density are generated.
*/
  F_IN.open( "jorge.inp", std::ifstream::in) ;

  while( getline( F_IN, line)){
    v( i) = std::stod( line) ;
    u( i) = std::sqrt( d1 - v(i)*v(i)) ;
    i++ ;
    }

  std::cout << " Read in correct values " << std::endl ;

  v1 = v ;
  n1 = v1.dot( v1) - N ;

  if ( n1 > d0){
    s2 = d2 ;
  } else {
    s2 = d1/d2 ;
    }

  std::cout << " Searching for proper Scaling " << std::endl ;

  for( i = 0; i < nbas; i++){
    v2(i) = v(i)/std::sqrt( v(i)*v(i) + s2*s2*u(i)*u(i)) ;
    }

  n2 = v2.dot(v2) - N ;

  for ( j = 0; j < 20; j++) {
    s3 = ( s1*n2 - s2*n1)/( n2 - n1) ;

    s1 = s2 ;
    n1 = n2 ;
    v1 = v2 ;

    s2 = s3 ;
    for( i = 0; i < nbas; i++){
      v2(i) = v(i)/std::sqrt( v(i)*v(i) + s2*s2*u(i)*u(i)) ;
      }
    n2 = v2.dot(v2) - N ;
    if ( std::abs( n2) < 1.0e-6) { break ;}
    }

  v = v2 ;
  for( i = 0; i < nbas; i++){
    u(i) = std::sqrt( d1 - v(i)*v(i)) ;
    }

  n1 = v.dot(v) ;

/*
  Now we build our rho and kappa guesses
*/
  for( i = 0; i < nbas; i++){
    p( i, i) = static_cast<cd>( v( i)*v( i)) ;
    k( i, i) = static_cast<cd>( v( i)*u( i)) ;
    }

  for( i = 0; i < nbas; i++){
    v1(i) = u(i)*u(i) ;
    }

  t = v1.asDiagonal() ;
  n2 = std::sqrt( t.determinant()) ;
  norm = n2*n2 ;
  std::cout << " norm " << norm << std::endl ;

  F_IN.close() ;
  u.resize( 0) ;
  v.resize( 0) ; 
  v1.resize( 0) ;
  v2.resize( 0) ;

  return ;

} ;
