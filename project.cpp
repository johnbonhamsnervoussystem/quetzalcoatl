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
  proj_HFB( com) ;

  return ;

} ;

void proj_HFB( common& com){
/*
  Projected HFB for Real wavefunctions

    - This will set up all the necessary routines to do various projection.
    - Set up integration grids

    - Let's start with Number projection as defualt.  We will do repreated diagonalization 
      of the Number projected Hmailtonian

    - This is only Generalized right now

    - Matrices are complex

  H - Core hamiltonian
  W - HFB wavefunction
  x - integration points
  wt - integration weights
*/

  int nbas = com.nbas() ;
  int maxit = com.mxscfit() ;
  int nalp = com.nalp() ; 
  double nele = static_cast<double>( com.nele()) ; 
  double mu = d0 ; 
  cd lshift = static_cast<cd>(com.lvlshft()), norm, olap ;
  diis_control diis_opt ;
  Eigen::MatrixXcd h ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXcd p, k, pt, kt ;
  Eigen::MatrixXcd t, u, prr, krr, I ;
  nbodyint<Eigen::MatrixXcd>* X ;
  std::vector<tei>* r12int ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w ;
  time_dbg proj_HFB_time = time_dbg("proj_HFB") ;

  t.resize( 2*nbas, 2*nbas) ;
  prr.resize( nbas, nbas) ;
  krr.resize( nbas, nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  D.resize( nbas, nbas) ;
  p.resize( 2*nbas, 2*nbas) ;
  k.resize( 2*nbas, 2*nbas) ;
  pt.resize( nbas, nbas) ;
  kt.resize( nbas, nbas) ;
  u.resize( 4*nbas, 4*nbas) ;
  w.moc.resize( 2*nbas, 2*nbas) ;
  w.eig.resize( 2*nbas) ;

  p.setZero() ;
  k.setZero() ;
  pt.setZero() ;
  kt.setZero() ;
  u.setZero() ;
  I.setIdentity() ;

/*
  Set diis options for now.  Wrap into a constructor later
*/
  diis_opt.set_switch( 2) ;
  diis_opt.do_diis = false ;
  diis_opt.diis_switch = false ;
  diis_opt.ndiisv = 5 ;
  diis_opt.diistype = 2 ;
  diis_opt.ediff_thr = 5.0e-3 ;
  
/*
  Build rho and kappa
*/
  int blahhhh = 1 ;
  real_SlaDet( com, blahhhh) ;
/*
  Rather than use the previous energy of each iteration to determine
  when the diis should turn on, we will use the hartree-fock energy
  as an upper bound to ensure we have passed that solution.
*/
  load_wfn( w) ;
  diis_opt.cntl_ref = w.e_scf ;

//  jorge_guess( pt, kt, Nalp, norm) ;
  thermal_guess( nalp, nbas, pt, kt) ;

/*
  I should be able to use rrhfb here.
*/
//  initialize( 2, 3, com.hamil(), com, h, X, r12int, nbas) ;

  w.moc.resize( 2*nbas, 2*nbas) ;
  w.eig.resize( 2*nbas) ;

  lshift = z4 ;
/*
  Build the particle number integration grid
*/
  for ( int qz=22; qz < 23; qz+=2) {
    trapezoid* ngrid = new trapezoid( d0, d2*pi, qz) ;
    h.setZero() ;
    initialize( 2, 1, com.hamil(), com, h, X, r12int, nbas) ;
    t = h.block( 0, 0, nbas, nbas) ;
/*
    prr = pt ;
    krr = kt ;
    Xring_shiekh_rr( nbas, t, X, w.moc, prr, krr, mu, lshift, ngrid, nele, diis_opt, maxit) ;
*/

    prr = pt ;
    krr = kt ;
    ring_shiekh_rr( nbas, t, X, w.moc, prr, krr, mu, lshift, ngrid, nele, diis_opt, maxit) ;
/*
  Reset the DIIS options just in case
*/
    diis_opt.set_switch( 2) ;
    diis_opt.do_diis = true ;
    diis_opt.diis_switch = false ;
    diis_opt.ndiisv = 5 ;
    diis_opt.diistype = 2 ;
    diis_opt.ediff_thr = 5.0e-3 ;
  
    w.moc.resize( 4*nbas, 4*nbas) ;
    w.eig.resize( 4*nbas) ;
    h.setZero() ;
    p.block( 0, 0, nbas, nbas) = pt ;
    p.block( nbas, nbas, nbas, nbas) = pt ;
    k.block( 0, nbas, nbas, nbas) = kt ;
    k.block( nbas, 0, nbas, nbas) = -kt ;
    t = k.inverse() ;
    u.block( 0, 0, 2*nbas, 2*nbas) = -p*t.conjugate() ;
    u.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
    u.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
    u.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = t*p ;
    olap = pfaffian_H( u) ;
    norm = z1/olap ;
    initialize( 2, 3, com.hamil(), com, h, X, r12int, nbas) ;
    t = h.block( 0, 0, 2*nbas, 2*nbas) ;
    ring_shiekh_cg( nbas, t, X, w.moc, p, k, mu, lshift, ngrid, nele, norm, diis_opt, maxit) ;
    }

  proj_HFB_time.end() ;

  return ;

} ;

template < class matrix>
void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  Real Restricted NHFB implementation using the Ring and Shiekh equations published in

  Sheikh, Javid A., and Peter Ring. "Symmetry-projected Hartree-Fock-Bogoliubov equations." 
  Nuclear Physics A 665 1-2 (2000): 71-91"

  This version needs improvement but it is not a priority right now.
*/
  int iter = 1, iter_N = 0 ;
  int in, inmb = ngrid->ns() ;
  double b_ll, b_ul ;
  double pnum = nele/d2 ;
  cd npar = static_cast<cd>(nele) ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  diis<cd> fdiis( 2*nbas, diis_cntrl.ndiisv, diis_cntrl.diistype) ;
  cd N = z0, olap ;
  cd energy, pphase, nphase ;
  cd s_theta, intE_g = z0, intx_g, fnx, w_phi, x_phi ;
  cd etr, gtr, dtr, d_phi ;
  cd shift ;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::MatrixXcd R, R_save, eigvec ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd tmp, Y_tmp, I ;
  Eigen::MatrixXcd C_phi, G_phi, D_phi ;
  Eigen::MatrixXcd r_phi, k_phi, R_phi ;
  Eigen::MatrixXcd Y_phi ;
  Eigen::MatrixXcd kbar_phi, t ;
  Eigen::MatrixXcd mu, mu_n, Heff, H ;

  time_dbg ring_shiekh_rr_projection_time = time_dbg("ring and shiekh rr") ;

/*
  Eigen Vectors to accumulate quantities
*/

  R.resize( 2*nbas, 2*nbas) ;
  R_save.resize( 2*nbas, 2*nbas) ;
  eigvec.resize( 2*nbas, 2*nbas) ;
  p_rho.resize( nbas, nbas) ;
  p_kappa.resize( nbas, nbas) ;
  epsilonN.resize( nbas, nbas) ;
  gammaN.resize( nbas, nbas) ;
  lambdaN.resize( nbas, nbas) ;
  FockN.resize( nbas, nbas) ;
  DeltaN.resize( nbas, nbas) ;
  Y_tmp.resize( nbas, nbas) ;
  tmp.resize( nbas, nbas) ;
  I.resize( nbas, nbas) ;
  C_phi.resize( nbas, nbas) ;
  G_phi.resize( nbas, nbas) ;
  R_phi.resize( nbas, nbas) ;
  D_phi.resize( nbas, nbas) ;
  r_phi.resize( nbas, nbas) ;
  k_phi.resize( nbas, nbas) ;
  Y_phi.resize( nbas, nbas) ;
  kbar_phi.resize( nbas, nbas) ;
  t.resize( nbas, nbas) ;
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
    Prepare level shifting density
  */

  R.block( 0, 0, nbas, nbas) = rho ;
  R.block( 0, nbas, nbas, nbas) = kappa ;
  R.block( nbas, 0, nbas, nbas) = kappa.adjoint() ;
  R.block( nbas, nbas, nbas, nbas) = I - rho.conjugate() ;

  do {
    /*
      Loop over the number projection grid and accumulate quantities 
      Clear our accumulation matrices
    */
    p_rho = rho ;
    p_kappa = kappa ;

    epsilonN.setZero() ;
    gammaN.setZero() ;
    lambdaN.setZero() ;
    FockN.setZero() ;
    DeltaN.setZero() ;
    Y_tmp.setZero() ;

    intx_g = z0 ;
    intE_g = z0 ;

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
      t = tmp.inverse() ;
      C_phi = pphase*t ;
      r_phi = C_phi*rho ;
      k_phi = C_phi*kappa ;
      kbar_phi = pphase*kappa*C_phi.conjugate() ;
      olap =  std::exp( zi*z2*fnx*static_cast<cd>(nbas))/C_phi.determinant() ;
      x_phi = d_phi*olap ;
      Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;
      intx_g += w_phi*x_phi ;
      Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

      /*
        Get the energy expression.
        Start by putting into the non orthogonal basis
      */

      W->contract( r_phi, k_phi) ;
      eigvec = W->getG() ;
      G_phi = eigvec.block( 0, 0, nbas, nbas) ;
      D_phi = eigvec.block( 0, nbas, nbas, nbas) ;

      t = h*r_phi ;
      etr = z2*t.trace() ;
      epsilonN = w_phi*x_phi*(Y_phi*etr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi) ;
      t = G_phi*r_phi ;
      gtr = t.trace() ;
      gammaN = w_phi*x_phi*(Y_phi*gtr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*G_phi*C_phi) ;
      t = D_phi*kbar_phi.adjoint() ;
      dtr = -t.trace()/z2 ;
      t = D_phi.transpose()*kbar_phi.conjugate() ;
      dtr -= t.trace()/z2 ;
      lambdaN = -w_phi*x_phi*(Y_phi*dtr + z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint()) ;

      FockN += epsilonN + gammaN + lambdaN ;
      DeltaN += w_phi*nphase*x_phi*C_phi*D_phi ;

/*
  Accumulate the overlap and projected energy
*/
      intE_g += w_phi*x_phi*( etr + gtr - dtr) ;
      }

    ngrid->set_s() ;

/*
  Now that everything has been accumulated, there are a few post loop modifications to make.
*/

    intE_g /= intx_g ;
    FockN /= intx_g ;
    DeltaN /= intx_g ;
    Y_tmp /= intx_g ;
    FockN -= intE_g*Y_tmp ;

    cnv_energy.push_back( intE_g) ;

    t = (FockN + FockN.adjoint())/z2 ;
    Heff.block( 0, 0, nbas, nbas) = t ;

    t = (DeltaN + DeltaN.transpose())/z2 ;
    Heff.block( 0, nbas, nbas, nbas) = t ;

    Heff.block( nbas, nbas, nbas, nbas) = -Heff.block( 0, 0, nbas, nbas).conjugate() ;
    Heff.block( nbas, 0, nbas, nbas) = Heff.block( 0, nbas, nbas, nbas).adjoint() ;

    if ( iter < 5 ){
      shift = lshift/static_cast<cd>(iter) ;
      }

    Heff += -shift*R ;

    Heff += shift*(mu_n - R) ;

    if ( diis_cntrl.do_diis){

      R_save = R ;

      }

  /*
    Adjust the chemical potential until the density gives the proper number
    of particles.
  */

  chemical_potential ( pnum, nbas, lambda, H_diag, R, mu, eigvec, Heff, rho, H) ;

  if ( diis_cntrl.do_diis){
    diis_cntrl.toggle( std::real(intE_g)) ;
    fdiis.update( R_save, H, diis_cntrl.diis_switch) ;
    if ( diis_cntrl.diis_switch ) {
      Heff = H ;
      chemical_potential ( pnum, nbas, lambda, H_diag, R, mu, eigvec, Heff, rho, H) ;
      }
    }

   kappa = R.block( 0, nbas, nbas, nbas) ;

/*
  Check for convergence
*/

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;

   if ( energy.real() < 1.0e-7) {
     std::cout << " Converged in the density after " << iter << " iterations. "  << std::endl ;
     std::cout << " Chemical potential: " << lambda << std::endl ;
     std::cout << " Energy " << intE_g << std::endl ;
     break ;
     }
 
 } while ( ++iter < maxit) ;

  ngrid->set_s() ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }

  ring_shiekh_rr_projection_time.end() ;

  return ;

} ;

template void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<Eigen::MatrixXcd>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control&, int& maxit) ;

template < class matrix>
void ring_shiekh_rr_K( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> V, Eigen::Ref<Eigen::MatrixXcd> U, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  Real Restricted KNHFB implementation.
*/

  time_dbg ring_shiekh_rr_K_projection_time = time_dbg("ring and shiekh rr K") ;

  ring_shiekh_rr_K_projection_time.end() ;

  return ;

} ;

template void ring_shiekh_rr_K( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd, trapezoid*&, double&, diis_control&, int&) ;

template < class matrix>
void Xring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  This is a working version of the routine however inefficient.  DO NOT ALTER
*/
  int iter = 1, iter_N = 0 ;
  int in, inmb = ngrid->ns() ;
  double b_ll, b_ul ;
  double pnum = nele/d2 ;
  cd npar = static_cast<cd>(nele) ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  diis<cd> fdiis( 2*nbas, diis_cntrl.ndiisv, diis_cntrl.diistype) ;
  cd N = z0, olap ;
  cd energy, pphase, nphase ;
  cd s_theta, intE_g = z0, intx_g, fnx, w_phi, x_phi, y_phi ;
  cd etr, gtr, dtr, d_phi ;
  cd shift, prev_energy ;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::MatrixXcd R, R_save, eigvec ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd /*DbaN,*/ tmp, Y_tmp, I ;
  Eigen::MatrixXcd C_phi, G_phi, D_phi ;
  Eigen::MatrixXcd r_phi, k_phi, R_phi ;
  Eigen::MatrixXcd Y_phi ;
  Eigen::MatrixXcd kbar_phi, t ;
  Eigen::MatrixXcd mu, mu_n, Heff, H ;

  time_dbg ring_shiekh_rr_projection_time = time_dbg("Xring and shiekh rr") ;

/*
  Eigen Vectors to accumulate quantities
*/

  R.resize( 2*nbas, 2*nbas) ;
  R_save.resize( 2*nbas, 2*nbas) ;
  eigvec.resize( 2*nbas, 2*nbas) ;
  p_rho.resize( nbas, nbas) ;
  p_kappa.resize( nbas, nbas) ;
  epsilonN.resize( nbas, nbas) ;
  gammaN.resize( nbas, nbas) ;
  lambdaN.resize( nbas, nbas) ;
  FockN.resize( nbas, nbas) ;
  DeltaN.resize( nbas, nbas) ;
  Y_tmp.resize( nbas, nbas) ;
  tmp.resize( nbas, nbas) ;
  I.resize( nbas, nbas) ;
  C_phi.resize( nbas, nbas) ;
  G_phi.resize( nbas, nbas) ;
  R_phi.resize( nbas, nbas) ;
  D_phi.resize( nbas, nbas) ;
  r_phi.resize( nbas, nbas) ;
  k_phi.resize( nbas, nbas) ;
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
    Prepare level shifting density
  */

  R.block( 0, 0, nbas, nbas) = rho ;
  R.block( 0, nbas, nbas, nbas) = kappa ;
  R.block( nbas, 0, nbas, nbas) = kappa.adjoint() ;
  R.block( nbas, nbas, nbas, nbas) = I - rho.conjugate() ;

  do {
    /*
      Loop over the number projection grid and accumulate quantities 
      Clear our accumulation matrices
    */
    p_rho = rho ;
    p_kappa = kappa ;

    epsilonN.setZero() ;
    gammaN.setZero() ;
    lambdaN.setZero() ;
    FockN.setZero() ;
    DeltaN.setZero() ;
    Y_tmp.setZero() ;

    prev_energy = intE_g ;
    intx_g = z0 ;
    intE_g = z0 ;
/*
  I don't want to do two loops because I think the same thing can be accomplished in one.
  However in the spirit of debugging, I will do two loops to match the equations Ring and 
  Shiekh wrote out as closely as possible.
*/

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
//      std::cout << " angle " << fnx << std::endl ;
      w_phi = ngrid->w() ;
//      std::cout << " weight " << w_phi << std::endl ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>(nele))/( z2*zpi) ;
      pphase = std::exp( z2*zi*fnx) ;
//      std::cout << " pphase " << pphase << std::endl ;

      tmp = I + rho*( pphase - z1) ;
//      print_mat( tmp, " I + p*( pphase  - z1 ) ") ;
      t = tmp.inverse() ;
//      print_mat( t, " tmp.inverse ") ;
      C_phi = pphase*t ;
//      print_mat( C_phi, " C (phi) ") ;
//      std::cout << " det(tmp.inverse()) " << t.determinant() << std::endl ;
//      std::cout << " fac1(nbas) " << std::exp( zi*fnx*static_cast<cd>(nbas)) << std::endl ;
      olap =  std::exp( zi*z2*fnx*static_cast<cd>(nbas))/C_phi.determinant() ;
//      std::cout << " olap " << olap << std::endl ;
      x_phi = d_phi*olap ;
//      std::cout << " x_phi CAJ " << std::exp( -zi*fnx*static_cast<cd>(nele))/(z2*zpi*t.determinant()) << std::endl ;
//      std::cout << " x_phi " << x_phi << std::endl ;
//      std::cout << " zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi " << std::endl << zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi << std::endl ;
      Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

      intx_g += w_phi*x_phi ;

   //   std::cout << std::endl << std::endl ;
      }

    ngrid->set_s() ;

    Y_tmp /= intx_g ;

//    std::cout << " intx_g " << intx_g << std::endl ;

//    print_mat( Y_tmp, "Y_tmp") ;

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
      t = tmp.inverse() ;
      C_phi = pphase*t ;
//      print_mat( C_phi, "C(phi)") ;
      r_phi = C_phi*rho ;
//      print_mat( r_phi, "p(phi)") ;
      k_phi = C_phi*kappa ;
//      print_mat( k_phi, "k(phi)") ;
      kbar_phi = pphase*kappa*C_phi.conjugate() ;
//      print_mat( kbar_phi, "kbar(phi)") ;
      olap =  std::exp( zi*z2*fnx*static_cast<cd>(nbas))/C_phi.determinant() ;
      x_phi = d_phi*olap ;
//      std::cout << " x_phi " << x_phi << std::endl ;
      y_phi = x_phi/intx_g ;
//      std::cout << " y_phi " << y_phi << std::endl ;
      Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi - Y_tmp ;

      /*
        Get the energy expression.
        Start by putting into the non orthogonal basis
      */

      W->contract( r_phi, k_phi) ;
      eigvec = W->getG() ;
      G_phi = eigvec.block( 0, 0, nbas, nbas) ;
//      print_mat( G_phi, " G(phi) ") ;
      D_phi = eigvec.block( 0, nbas, nbas, nbas) ;
//      print_mat( D_phi, " D(phi) ") ;

      t = h*r_phi ;
      etr = z2*t.trace() ;
//      std::cout << " 2*tr(hp) " << etr << std::endl ;
//      t = y_phi*Y_phi*etr ;
//      print_mat( t, " y_phi*Y_phi*etr ") ;
//      t = y_phi*( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi ;
//      print_mat( t, " y_phi*( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi ") ;
      epsilonN = w_phi*y_phi*(Y_phi*etr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*h*C_phi) ;
      t = G_phi*r_phi ;
      gtr = t.trace() ;
      gammaN = w_phi*y_phi*(Y_phi*gtr + ( I - z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*r_phi)*G_phi*C_phi) ;
      t = D_phi*kbar_phi.adjoint() ;
      dtr = -t.trace()/z2 ;
      t = D_phi.transpose()*kbar_phi.conjugate() ;
      dtr -= t.trace()/z2 ;
//      std::cout << " -(tr(D^{T}*kbar_phi^{*}) + tr(D*kbar_phi^{t}))/2 " << dtr << std::endl ;
//      std::cout << " y_phi " << y_phi << std::endl ;
//      print_mat( Y_phi, "Y(phi)") ;
//      t = y_phi*Y_phi*dtr ;
//      print_mat( t, " y_phi*Y_phi*dtr ") ;
//      t = y_phi*z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint() ;
//      print_mat( t, " y_phi*z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint() ") ;
      lambdaN = -w_phi*y_phi*(Y_phi*dtr + z2*zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi*D_phi*kbar_phi.adjoint()) ;

      FockN += epsilonN + gammaN + lambdaN ;
//      t = nphase*y_phi*C_phi*D_phi ;
//      print_mat( t, " nphase*y_phi*C_phi*D_phi ") ;
      DeltaN += w_phi*nphase*y_phi*C_phi*D_phi ;
//      DbaN = w_phi*nphase*y_phi*C_phi*D_phi.transpose() ;

/*
  Accumulate the overlap and projected energy
*/
      intE_g += w_phi*y_phi*( etr + gtr - dtr) ;
      }

    ngrid->set_s() ;

/*
    t = (epsilonN + epsilonN.adjoint())/z2 ;
    print_mat( t, " FN") ;

    t = (gammaN + gammaN.adjoint())/z2 ;
    print_mat( t, " GN") ;

    t = (lambdaN + lambdaN.adjoint())/z2 ;
    print_mat( t, " LN") ;

    t = (DeltaN - DbaN.transpose())/z2 ;
    print_mat( t, " DabN") ;

    t = (DeltaN + DeltaN.transpose())/z2 ;
    print_mat( t, " DabN") ;
*/

    cnv_energy.push_back( intE_g) ;

    t = (FockN + FockN.adjoint())/z2 ;
//    print_mat( t, " FockN ") ;
    Heff.block( 0, 0, nbas, nbas) = t ;

    t = (DeltaN + DeltaN.transpose())/z2 ;
    Heff.block( 0, nbas, nbas, nbas) = t ;

    Heff.block( nbas, nbas, nbas, nbas) = -Heff.block( 0, 0, nbas, nbas).conjugate() ;
    Heff.block( nbas, 0, nbas, nbas) = Heff.block( 0, nbas, nbas, nbas).adjoint() ;

    if ( iter < 5 ){
      shift = lshift/static_cast<cd>(iter) ;
      }

    Heff += -shift*R ;

    Heff += shift*(mu_n - R) ;

    if ( diis_cntrl.do_diis){

      R_save = R ;

      }

  /*
    Adjust the chemical potential until the density gives the proper number
    of particles.
  */

    if ( true ){
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
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        } while ( ++iter_N < 30 && static_cast<double>(N.real()) > pnum) ;
  
      iter_N = 0 ;
  
      do {
        b_ul += d3 ;
        H = Heff + static_cast<cd>(b_ul)*mu ;
        H_diag.compute( H) ;
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        } while ( ++iter_N < 30 && static_cast<double>(N.real()) < pnum) ;
  
      iter_N = 0 ;
  
      lambda = ( b_ll + b_ul)/d2 ;
  
      while ( iter_N++ < 30) {
        H = Heff + static_cast<cd>(lambda)*mu ;
        H_diag.compute( H) ;
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        if ( std::abs(static_cast<double>(N.real()) - pnum) < 1.0e-5){
          break ;
        } else {
  
/*
  Bisection method
*/

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
        H_diag.compute( Heff) ;
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        }

    if ( diis_cntrl.do_diis){
      diis_cntrl.toggle( std::real(intE_g)) ;
      fdiis.update( R_save, H, diis_cntrl.diis_switch) ;
      }

   kappa = R.block( 0, nbas, nbas, nbas) ;

/*
  Check for convergence
*/

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;

   if ( energy.real() < 1.0e-7) {
     std::cout << " Converged in the density after " << iter << " iterations. "  << std::endl ;
     std::cout << " Chemical potential: " << lambda << std::endl ;
     std::cout << " Energy " << intE_g << std::endl ;
     break ;
     }
 
 } while ( ++iter < maxit) ;

  ngrid->set_s() ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }

  ring_shiekh_rr_projection_time.end() ;

  return ;

} ;

template void Xring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<Eigen::MatrixXcd>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control&, int& maxit) ;

template < class matrix>
void ring_shiekh_cg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, cd &nrm, diis_control& diis_cntrl, int& maxit) {

/*
  Various debugging for the number projected HFB.
*/

  int iter = 1, iter_N = 0, dbas = 2*nbas ;
  int in, inmb = ngrid->ns() ;
  double b_ll, b_ul ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  diis<cd> fdiis( 4*nbas, diis_cntrl.ndiisv, diis_cntrl.diistype) ;
  cd N = z0, olap, shift ;
  cd energy, pphase, nphase, d_phi ;
  cd s_theta, intE_g, intx_g, fnx, w_phi, x_phi ;
  cd etr, gtr, dtr, cds_scr ;
  std::vector<cd> cnv_den, cnv_energy ;
  Eigen::MatrixXcd R, R_save, eigvec ;
  Eigen::MatrixXcd p_rho, p_kappa ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd tmp, Y_tmp, I ;
  Eigen::MatrixXcd C_phi, G_phi, D_phi ;
  Eigen::MatrixXcd r_phi, k_phi, R_phi ;
  Eigen::MatrixXcd X_phi, Y_phi ;
  Eigen::MatrixXcd kbar_phi, t, t1, kappa_i ;
  Eigen::MatrixXcd mu, mu_n, Heff, H ;

  time_dbg ring_shiekh_projection_time = time_dbg("ring and shiekh cg") ;

/*
  Eigen Vectors to accumulate quantities
*/

  R.resize( 4*nbas, 4*nbas) ;
  R_save.resize( 4*nbas, 4*nbas) ;
  eigvec.resize( 4*nbas, 4*nbas) ;
  p_rho.resize( dbas, dbas) ;
  p_kappa.resize( dbas, dbas) ;
  epsilonN.resize( dbas, dbas) ;
  gammaN.resize( dbas, dbas) ;
  lambdaN.resize( dbas, dbas) ;
  FockN.resize( dbas, dbas) ;
  DeltaN.resize( dbas, dbas) ;
  Y_tmp.resize( dbas, dbas) ;
  tmp.resize( dbas, dbas) ;
  I.resize( dbas, dbas) ;
  C_phi.resize( dbas, dbas) ;
  G_phi.resize( dbas, dbas) ;
  R_phi.resize( dbas, dbas) ;
  D_phi.resize( dbas, dbas) ;
  r_phi.resize( dbas, dbas) ;
  k_phi.resize( dbas, dbas) ;
  kappa_i.resize( dbas, dbas) ;
  X_phi.resize( dbas, dbas) ;
  Y_phi.resize( dbas, dbas) ;
  kbar_phi.resize( dbas, dbas) ;
  t.resize( dbas, dbas) ;
  t1.resize( dbas, dbas) ;
  mu.resize( 4*nbas, 4*nbas) ;
  mu_n.resize( 4*nbas, 4*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;
  H.resize( 4*nbas, 4*nbas) ;

  I.setIdentity() ;
  mu.setZero() ;
  mu.block( 0, 0, dbas, dbas) = -I ;
  mu.block( dbas, dbas, dbas, dbas) = I ;
  mu_n.setIdentity() ;

  /*
    Initialize the generalized density
  */

  R.block( 0, 0, dbas, dbas) = rho ;
  R.block( 0, dbas, dbas, dbas) = kappa ;
  R.block( dbas, 0, dbas, dbas) = kappa.adjoint() ;
  R.block( dbas, dbas, dbas, dbas) = I - rho.conjugate() ;

  do {
    /*
      Loop over the number projection grid and accumulate quantities 
      Clear our accumulation matrices
    */
    p_rho = rho ;
    p_kappa = kappa ;
    kappa_i = kappa.inverse() ;

    epsilonN.setZero() ;
    gammaN.setZero() ;
    lambdaN.setZero() ;
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
      nphase = std::exp( -z2*zi*fnx) ;

      /*
        Build our rotated quantities
      */

      R_phi = I*std::exp( zi*fnx) ;
      tmp = I + rho*( pphase - z1) ;
      t = tmp.inverse() ;
      C_phi = pphase*t ;
      r_phi = C_phi*rho ;
      k_phi = C_phi*kappa ;
      kbar_phi = pphase*kappa*C_phi.conjugate() ;
      t = R_phi*kappa ;
      t1 = t.inverse() ;
      eigvec.block( 0, 0, dbas, dbas) = -R_phi*rho*t1.conjugate() ;
      eigvec.block( 0, dbas, dbas, dbas) = -I ;
      eigvec.block( dbas, 0, dbas, dbas) = I ;
      eigvec.block( dbas, dbas, dbas, dbas) = kappa_i*rho ;
      olap = nrm*pfaffian_H( eigvec) ;
      x_phi = d_phi*olap ;

      Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

      intx_g += w_phi*x_phi ;
      Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

      W->contract( r_phi, k_phi) ;
      eigvec = W->getG() ;
      G_phi = eigvec.block( 0, 0, dbas, dbas) ;
      D_phi = eigvec.block( 0, dbas, dbas, dbas) ;

      cds_scr = z2*zi*std::exp( -zi*fnx)*std::sin( fnx) ;
      t = h*r_phi ;
      etr = t.trace() ;
      epsilonN = w_phi*x_phi*(Y_phi*etr + ( I - cds_scr*r_phi)*h*C_phi) ;
      t = G_phi*r_phi ;
      gtr = t.trace() ;
      gammaN = w_phi*x_phi*(Y_phi*gtr/z2 + ( I - cds_scr*r_phi)*G_phi*C_phi) ;
      t = D_phi*kbar_phi.adjoint() ;
      dtr = -t.trace()/z2 ;
      t = D_phi.transpose()*kbar_phi.conjugate() ;
      dtr -= t.trace()/z2 ;
      lambdaN = -w_phi*x_phi*(Y_phi*dtr/z2 + cds_scr*C_phi*D_phi*kbar_phi.adjoint()) ;

      FockN += epsilonN + gammaN + lambdaN ;
      DeltaN += w_phi*nphase*x_phi*C_phi*D_phi ;

/*
  Accumulate the overlap and projected energy
*/
      intE_g += w_phi*x_phi*( etr + gtr/z2 - dtr/z2) ;
      }

    ngrid->set_s() ;

/*
  Post Processing of the accumulated quantities
*/

    intE_g /= intx_g ;
    FockN /= intx_g ;
    DeltaN /= intx_g ;
    Y_tmp /= intx_g ;
    FockN -= intE_g*Y_tmp ;

    cnv_energy.push_back( intE_g) ;

    t = (FockN + FockN.adjoint())/z2 ;
    Heff.block( 0, 0, dbas, dbas) = t ;
    t = (DeltaN - DeltaN.transpose())/z2 ;
    Heff.block( 0, dbas, dbas, dbas) = t ;

    Heff.block( dbas, dbas, dbas, dbas) = -Heff.block( 0, 0, dbas, dbas).conjugate() ;
    Heff.block( dbas, 0, dbas, dbas) = -Heff.block( 0, dbas, dbas, dbas).conjugate() ;

    if ( iter < 5){
      shift = lshift/static_cast<cd>(iter) ;
      }

    Heff += -shift*R ;

    Heff += shift*(mu_n - R) ;
/*
  Adjust the chemical potential until the density gives the proper number
  of particles.
*/

    if ( diis_cntrl.do_diis){

      R_save = R ;

      }

   chemical_potential ( nele, dbas, lambda, H_diag, R, mu, eigvec, Heff, rho, H) ;

   if ( diis_cntrl.do_diis){
     diis_cntrl.toggle( std::real(intE_g)) ;
     fdiis.update( R_save, H, diis_cntrl.diis_switch) ;
      if ( false ) {
/*      if ( diis_cntrl.diis_switch ) {
  For some reason this second adjustment of the chemical potential leads to 
  convergence to a higher energy solution.  Avoid it for now.
*/
        Heff = H ;
        chemical_potential ( nele, dbas, lambda, H_diag, R, mu, eigvec, Heff, rho, H) ;
        }
     }

   kappa = R.block( 0, dbas, dbas, dbas) ;
   t = eigvec.block( dbas, 0, dbas, dbas).conjugate()*eigvec.block( dbas, 0, dbas, dbas).transpose() ;
   nrm = std::sqrt( t.determinant()) ;

/*
  Check for convergence
*/

   t = rho - p_rho ;
   energy = t.norm() ;
   t = kappa - p_kappa ;
   energy += t.norm() ;
   cnv_den.push_back( energy) ;

   if ( energy.real() < 1.0e-7 ) {
     std::cout << " Converged in the density after " << iter << " iterations." << std::endl ;
     std::cout << " Chemical potential: " << lambda << std::endl ;
     std::cout << " Energy " << intE_g << std::endl ;
     break ;
     }
 
 } while ( ++iter < maxit) ;


  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }


  ring_shiekh_projection_time.end() ;

  return ;

} ;

template void ring_shiekh_cg( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd, trapezoid*&, double&, cd&, diis_control&, int&) ;

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

template < class matrix>
void projected_wavefunction_energy( const int& nele, const int& nbas, nbodyint<matrix>* W, trapezoid*& ngrid, const Eigen::Ref<Eigen::MatrixXcd> h, Eigen::MatrixXcd rho, Eigen::MatrixXcd kappa, const Eigen::Ref<Eigen::MatrixXcd> kappa_i, Eigen::Ref<Eigen::MatrixXcd> R_phi, Eigen::Ref<Eigen::MatrixXcd> r_phi, Eigen::Ref<Eigen::MatrixXcd> k_phi, Eigen::Ref<Eigen::MatrixXcd> kbar_phi, Eigen::Ref<Eigen::MatrixXcd> I, Eigen::Ref<Eigen::MatrixXcd> scr1, Eigen::Ref<Eigen::MatrixXcd> scr2, Eigen::Ref<Eigen::MatrixXcd> scr3, cd& intx_g, cd& intE_g) {
/*
  This is a wrapper for convienence and readability and so most everything is passed preallocated.
  Given a density, pairing matrix, two electron interaction, and integration grids, calculate the energy and overlap.
*/

  int in, inmb = ngrid->ns()  ;
  cd fnx, w_phi, x_phi, d_phi, pphase, nphase, olap, nrm, etr ;

  for ( in = 0; in < inmb; in++){
    fnx = static_cast<cd>( ngrid->q()) ;
    w_phi = ngrid->w() ;
    d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
    pphase = std::exp( z2*zi*fnx) ;
    nphase = std::exp( -z2*zi*fnx) ;

    /*
      Build our rotated quantities
    */

    R_phi = I*std::exp( zi*fnx) ;
    scr1 = I + rho*( pphase - z1) ;
    scr2 = scr1.inverse() ;
    scr1 = pphase*scr2 ;
    r_phi = scr1*rho ;
    k_phi = scr1*kappa ;
    kbar_phi = pphase*kappa*scr1.conjugate() ;
    scr1 = R_phi*kappa ;
    scr2 = scr1.inverse() ;
    scr3.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*scr2.conjugate() ;
    scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
    scr3.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;
    scr3.block(  2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;

    if ( in == 0){
      /* This is our unrotated wavefunction and is the norm */
      olap = pfaffian_H( scr3) ;
      nrm = z1/olap ;
      olap = z1 ;
    } else {
      olap = nrm*pfaffian_H( scr3) ;
      }

    x_phi = d_phi*olap ;

    intx_g += w_phi*x_phi ;

    W->contract( r_phi, k_phi) ;
    scr3 = W->getG() ;

    scr1 = h*r_phi ;
    etr = scr1.trace() ;
    scr1 = scr3.block( 0, 0, 2*nbas, 2*nbas) ;
    scr2 = scr1*r_phi ;
    etr += scr2.trace() ;
    scr1 = scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
    scr2 = scr1*kbar_phi.adjoint() ;
    etr -= scr2.trace()/z2 ;
    scr2 = scr1.transpose()*kbar_phi.conjugate() ;
    etr -= scr2.trace()/z2 ;

/*
  Accumulate the overlap and projected energy
*/
    intE_g += w_phi*x_phi*etr ;
    }

  ngrid->set_s() ;

  return ;

} ;

template void projected_wavefunction_energy( const int&, const int&, nbodyint<Eigen::MatrixXcd>*, trapezoid*&, const Eigen::Ref<Eigen::MatrixXcd>, Eigen::MatrixXcd, Eigen::MatrixXcd, const Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, cd&, cd&) ;

void chemical_potential ( const int& nele, const int& nbas, double& lambda, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag, Eigen::Ref<Eigen::MatrixXcd> R, Eigen::Ref<Eigen::MatrixXcd> mu, Eigen::Ref<Eigen::MatrixXcd> eigvec, const Eigen::Ref<Eigen::MatrixXcd> Heff, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> H) {
/*
  This is a wrapper for the bisection method for locating chemical potential.
  This will make the code more readable.  Consequently, nearly everything is passed as
  a pre-allocated matrix so that memory is not allocated and deallocated over and over.
*/

      int iter_N = 0 ;

      double b_ul = lambda ;
      double b_ll = lambda ;
      cd N ;

      /*
         Set some initial limits
      */

      do {
        b_ll -= d3 ;
        H = Heff + static_cast<cd>(b_ll)*mu ;
        H_diag.compute( H) ;
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        } while ( ++iter_N < 30 && static_cast<double>(N.real()) > static_cast<double>(nele)) ;
    
      iter_N = 0 ;
    
      do {
        b_ul += d3 ;
        H = Heff + static_cast<cd>(b_ul)*mu ;
        H_diag.compute( H) ;
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        } while ( ++iter_N < 30 && static_cast<double>(N.real()) < static_cast<double>(nele)) ;
    
      iter_N = 0 ;
    
      lambda = ( b_ll + b_ul)/d2 ;
    
      while ( iter_N++ < 30) {
        H = Heff + static_cast<cd>(lambda)*mu ;
        H_diag.compute( H) ;
        eigvec = H_diag.eigenvectors() ;
        R = eigvec.block( 0, 0, 2*nbas, nbas)*eigvec.block( 0, 0, 2*nbas, nbas).adjoint() ;
        rho = R.block( 0, 0, nbas, nbas) ;
        N = rho.trace() ;
        if ( std::abs(static_cast<double>(N.real()) - static_cast<double>(nele)) < 1.0e-5){
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
        }

  return ;

} ;

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

  v1 = v ;
  n1 = v1.dot( v1) - N ;

  if ( n1 > d0){
    s2 = d2 ;
  } else {
    s2 = d1/d2 ;
    }

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
