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
#include <string>
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
  Eigen::MatrixXd A, B ;
  Eigen::MatrixXcd h, U, V ;
  Eigen::MatrixXcd p, k, pt, kt ;
  Eigen::MatrixXcd t, u, prr, krr, I ;
  nbodyint<Eigen::MatrixXcd>* X ;
  std::vector<tei>* r12int ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w ;
  time_dbg proj_HFB_time = time_dbg("proj_HFB") ;

  A.resize( nbas, nbas) ;
  B.resize( nbas, nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  prr.resize( nbas, nbas) ;
  krr.resize( nbas, nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  p.resize( 2*nbas, 2*nbas) ;
  k.resize( 2*nbas, 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;
  V.resize( 2*nbas, 2*nbas) ;
  pt.resize( nbas, nbas) ;
  kt.resize( nbas, nbas) ;
  u.resize( 2*nbas, nbas) ;
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
  diis_opt.ediff_thr = 5.0e-2 ;
  
/*
  Build rho and kappa
  int blahhhh = 1 ;
  real_SlaDet( com, blahhhh) ;
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
  for ( int qz=2; qz < 3; qz+=2) {
    trapezoid* ngrid = new trapezoid( d0, d2*pi, qz) ;
/*
    h.setZero() ;
    initialize( 2, 1, com.hamil(), com, h, X, r12int, nbas) ;
    t = h.block( 0, 0, nbas, nbas) ;

    prr = pt ;
    krr = kt ;
    ring_shiekh_rr( nbas, t, X, w.eig, w.moc, prr, krr, mu, lshift, ngrid, nele, diis_opt, maxit) ;

    ring_shiekh_rr( nbas, t, X, w.eig, w.moc, prr, krr, mu, lshift, ngrid, nele, diis_opt, maxit) ;
*/
    p.setZero() ;
    k.setZero() ;
    p.block( 0, 0, nbas, nbas) = pt ;
    p.block( nbas, nbas, nbas, nbas) = pt ;
    k.block( 0, nbas, nbas, nbas) = kt ;
    k.block( nbas, 0, nbas, nbas) = -kt ;
    initialize( 2, 3, com.hamil(), com, h, X, r12int, nbas) ;
    print_mat( h, " core hamiltonian ") ;
    print_mat( p, " density ") ;
    print_mat( k, " pairing ") ;
    mu = d0 ;
    general_derivative_testing( nbas, h, X, w.eig, w.moc, p, k, mu, lshift, ngrid, nele, diis_opt, maxit) ;

    return ;

/*
  Test the projection energy routine
*/
    mu = d0 ;
    lshift = z4 ;
/*
  Generate a random unitary matrix to mix the solution and search for a conjugate broken 
  solution

    U.setRandom() ;
    D = U.adjoint() - U ;

    ahm_exp( D, U, nbas, 0) ;
    U /= z10 ;
     U_{1} = U_{0} = V_{0}^{*}Z_{1}^{*}
    u.block( nbas, 0, nbas, nbas) = w.moc.block( nbas, 0, nbas, nbas).conjugate() + w.moc.block( 0, 0, nbas, nbas)*U ;
    u.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nbas).conjugate() + w.moc.block( nbas, 0, nbas, nbas)*U ;
    Reorthogonalize our vectors 
    prr = u.adjoint()*u ;
    print_mat( prr, " C^{dagger}C ") ;
    canort( prr, krr, nbas) ;
    U = krr.adjoint()*prr*krr ;
    print_mat( U, " U^{dagger}*C*U ") ;
    w.moc.block( 0, 0, 2*nbas, nbas) = u.block( 0, 0, 2*nbas, nbas)*krr ;
    U = w.moc.block( 0, 0, 2*nbas, nbas).adjoint()*w.moc.block( 0, 0, 2*nbas, nbas) ;
    print_mat( w.eig, " Energy Levels ") ;
    thermal_occm( nbas, nalp, w.eig, temp, A, B) ;
*/
//    homo_lumo_mix_c( w.moc, nbas, pi/d4) ;

    t = w.moc.block( 0, 0, 2*nbas, nbas).adjoint()*w.moc.block( 0, 0, 2*nbas, nbas) ;
    print_mat( t, " C^{t}C ") ;
    V.setZero() ;
    U.setZero() ;
    V.block( 0, 0, nbas, nbas) = w.moc.block( 0, 0, nbas, nbas) ;
    V.block( nbas, nbas, nbas, nbas) = w.moc.block( 0, 0, nbas, nbas) ;
    U.block( 0, nbas, nbas, nbas) = -w.moc.block( nbas, 0, nbas, nbas) ;
    U.block( nbas, 0, nbas, nbas) = w.moc.block( nbas, 0, nbas, nbas) ;

    t = V*V.adjoint() ;
    print_mat( t, " Density ") ;
    t = V*U.adjoint() ;
    print_mat( t, " pairing ") ;

    initialize( 2, 3, com.hamil(), com, h, X, r12int, nbas) ;
    print_mat( h, " Core Hamiltonian ") ;

    ring_shiekh_K( nbas, h, X, V, U, mu, lshift, ngrid, nele, diis_opt, maxit) ;
/*
  Reset the DIIS options just in case
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
*/
    }

  proj_HFB_time.end() ;

  return ;

} ;

template <class matrix>
void build_projected_hamiltonian_energy( const int& nele, const int& nbas, nbodyint<matrix>* W, trapezoid*& ngrid, const Eigen::Ref<Eigen::MatrixXcd> h, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, const Eigen::Ref<Eigen::MatrixXcd> kappa_i, Eigen::Ref<Eigen::MatrixXcd> R_phi, matrix& r_phi, matrix& k_phi, Eigen::Ref<Eigen::MatrixXcd> kbar_phi, Eigen::Ref<Eigen::MatrixXcd> C_phi, Eigen::Ref<Eigen::MatrixXcd> G_phi, Eigen::Ref<Eigen::MatrixXcd> D_phi, Eigen::Ref<Eigen::MatrixXcd> I, Eigen::Ref<Eigen::MatrixXcd> epsilonN, Eigen::Ref<Eigen::MatrixXcd> gammaN, Eigen::Ref<Eigen::MatrixXcd> lambdaN, Eigen::Ref<Eigen::MatrixXcd> FockN, Eigen::Ref<Eigen::MatrixXcd> DeltaN, Eigen::Ref<Eigen::MatrixXcd> Y_tmp, Eigen::Ref<Eigen::MatrixXcd> Y_phi, Eigen::Ref<Eigen::MatrixXcd> Heff, Eigen::Ref<Eigen::MatrixXcd> t, Eigen::Ref<Eigen::MatrixXcd> tmp, cd& intx_g, cd& intE_g) {
/*
  This is a wrapper for convienence and readability and so most everything is passed preallocated.
  Given a density, pairing matrix, two electron interaction, and integration grids, calculate the energy and overlap
  and generate the effective-Hamiltonian for projected methods.

  Super wasteful right now for debugging.  Will clean it later.
*/

  int in, dbas = 2*nbas, inmb = ngrid->ns() ;
  cd fnx, w_phi, x_phi, d_phi, pphase, nphase, olap, nrm ;
  cd cds_scr, etr, dtr, gtr ;

  time_dbg projected_wavefunction_energy_time = time_dbg("projected wavefunction energy") ;

  intx_g = z0 ;
  intE_g = z0 ;
  nrm = z0 ;

  ngrid->set_s() ;

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
     tmp = t.inverse() ;
     Heff.block( 0, 0, dbas, dbas) = -R_phi*rho*tmp.conjugate() ;
     Heff.block( 0, dbas, dbas, dbas) = -I ;
     Heff.block( dbas, 0, dbas, dbas) = I ;
     Heff.block( dbas, dbas, dbas, dbas) = kappa_i*rho ;
     if ( in == 0 ){
       nrm = z1/pfaffian_H( Heff) ;
       olap = nrm/nrm ;
     } else {
       olap = nrm*pfaffian_H( Heff) ;
       }
     x_phi = d_phi*olap ;

     Y_tmp += w_phi*zi*x_phi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

     intx_g += w_phi*x_phi ;
     Y_phi = zi*std::exp( -zi*fnx)*std::sin( fnx)*C_phi ;

     W->contract( r_phi, k_phi) ;
     Heff = W->getG() ;
     G_phi = Heff.block( 0, 0, dbas, dbas) ;
     D_phi = Heff.block( 0, dbas, dbas, dbas) ;

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

    FockN /= intx_g ;
    DeltaN /= intx_g ;
    Y_tmp /= intx_g ;
    FockN -= intE_g*Y_tmp/intx_g ;

    t = (FockN + FockN.adjoint())/z2 ;
    Heff.block( 0, 0, dbas, dbas) = t ;
    t = (DeltaN - DeltaN.transpose())/z2 ;
    Heff.block( 0, dbas, dbas, dbas) = t ;

    Heff.block( dbas, dbas, dbas, dbas) = -Heff.block( 0, 0, dbas, dbas).conjugate() ;
    Heff.block( dbas, 0, dbas, dbas) = -Heff.block( 0, dbas, dbas, dbas).conjugate() ;

  ngrid->set_s() ;

  projected_wavefunction_energy_time.end() ;

  return ;

} ;

template void build_projected_hamiltonian_energy( const int&, const int&, nbodyint<Eigen::MatrixXcd>*, trapezoid*&, const Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, const Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, cd&, cd&) ;


template < class matrix>
void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::VectorXd> eigval, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  Real Restricted NHFB implementation using the Ring and Shiekh equations published in

  Sheikh, Javid A., and Peter Ring. "Symmetry-projected Hartree-Fock-Bogoliubov equations." 
  Nuclear Physics A 665 1-2 (2000): 71-91"

  This version needs improvement but it is not a priority right now.
*/
  int iter = 1 ;
  int in, inmb = ngrid->ns() ;
  double pnum = nele/d2 ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
  diis<cd> fdiis( 2*nbas, diis_cntrl.ndiisv, diis_cntrl.diistype) ;
  cd olap ;
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

  /* Grab the last thing in H_diag */
  eigval = H_diag.eigenvalues().real() ;
  c = H_diag.eigenvectors() ;
/* 
  Print the final density and pairing for debugging
*/
  print_mat( rho, " Final Density ") ;
  print_mat( kappa, " Final Pairing ") ;

  for( unsigned int xq=0; xq < cnv_den.size(); xq++){
    std::cout << cnv_den[xq] << "  " << cnv_energy[xq] << std::endl ;
    }

  ring_shiekh_rr_projection_time.end() ;

  return ;

} ;

template void ring_shiekh_rr( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::VectorXd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd, trapezoid*&, double&, diis_control&, int&) ;

template < class matrix>
void ring_shiekh_K( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> V, Eigen::Ref<Eigen::MatrixXcd> U, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  KNHFB implementation.
*/

  int i ;
  cd PE, PO, CI_S ;
  Eigen::MatrixXcd rho, kappa, kappa_i ;
  Eigen::MatrixXcd R_phi, r_phi, k_phi ;
  Eigen::MatrixXcd kbar_phi, I, C_phi ;
  Eigen::MatrixXcd FockN, DeltaN, G_phi, D_phi ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, Heff ;
  Eigen::MatrixXcd t, tmp, Y_tmp, Y_phi ;
  Eigen::MatrixXcd Vx, Ux ;
  Eigen::VectorXcd e_vec( 2) ;
  Eigen::MatrixXcd HCI( 2, 2), SCI( 2, 2), CIVEC( 2, 2) ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> H_KCI ;
  std::vector<Eigen::MatrixXcd> h_vec ;

  rho.resize( 2*nbas, 2*nbas) ;
  kappa.resize( 2*nbas, 2*nbas) ;
  kappa_i.resize( 2*nbas, 2*nbas) ;
  Vx.resize( 2*nbas, 2*nbas) ;
  Ux.resize( 2*nbas, 2*nbas) ;
  R_phi.resize( 2*nbas, 2*nbas) ;
  r_phi.resize( 2*nbas, 2*nbas) ;
  k_phi.resize( 2*nbas, 2*nbas) ;
  C_phi.resize( 2*nbas, 2*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  tmp.resize( 2*nbas, 2*nbas) ;
  kbar_phi.resize( 2*nbas, 2*nbas) ;
  I.resize( 2*nbas, 2*nbas) ;
  epsilonN.resize( 2*nbas, 2*nbas) ;
  gammaN.resize( 2*nbas, 2*nbas) ;
  lambdaN.resize( 2*nbas, 2*nbas) ;
  FockN.resize( 2*nbas, 2*nbas) ;
  DeltaN.resize( 2*nbas, 2*nbas) ;
  Y_tmp.resize( 2*nbas, 2*nbas) ;
  Y_phi.resize( 2*nbas, 2*nbas) ;
  G_phi.resize( 2*nbas, 2*nbas) ;
  D_phi.resize( 2*nbas, 2*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;

  I.setIdentity() ;

  time_dbg ring_shiekh_K_projection_time = time_dbg("ring and shiekh K") ;

  for ( i = 0; i < 4; i++){
    if ( i == 0 ){
      std::cout << " <0|R|0> " << std::endl ;
      rho = V*V.adjoint() ;
      kappa = V*U.adjoint() ;
    } else if ( i == 2 ){
      std::cout << " <0|RK|0> " << std::endl ;
      Vx = V.conjugate() ;
      rho = Vx*V.adjoint() ;
      kappa = Vx*U.adjoint() ;
    } else if ( i == 1 ){
      std::cout << " <0|KR|0> " << std::endl ;
      Vx = V.conjugate() ;
      Ux = U.conjugate() ;
      rho = V*Vx.adjoint() ;
      kappa = V*Ux.adjoint() ;
    } else if ( i == 3 ){
      std::cout << " <0|KRK|0> " << std::endl ;
      Vx = V.conjugate() ;
      Ux = U.conjugate() ;
      rho = Vx*Vx.adjoint() ;
      kappa = Vx*Ux.adjoint() ;
      }
/*
    print_mat( rho, " rho ") ;
    print_mat( kappa, " kappa ") ;
*/
    R_phi.setZero() ;
    r_phi.setZero() ;
    k_phi.setZero() ;
    kbar_phi.setZero() ;
    C_phi.setZero() ;
    epsilonN.setZero() ;
    gammaN.setZero() ;
    lambdaN.setZero() ;
    FockN.setZero() ;
    DeltaN.setZero() ;
    Y_tmp.setZero() ;
    Y_phi.setZero() ;
    G_phi.setZero() ;
    D_phi.setZero() ;
    Heff.setZero() ;
    t.setZero() ;
    tmp.setZero() ;
    kappa_i = kappa.inverse() ;

//    projected_wavefunction_energy( nele, nbas, W, ngrid, h, rho, kappa, kappa_i, R_phi, r_phi, k_phi, kbar_phi, I, epsilonN, gammaN, Heff, PO, PE) ;
    build_projected_hamiltonian_energy( nele, nbas, W, ngrid, h, rho, kappa, kappa_i, R_phi, r_phi, k_phi, kbar_phi, C_phi, G_phi, D_phi, I, epsilonN, gammaN, lambdaN, FockN, DeltaN, Y_tmp, Y_phi, Heff, t, tmp, PO, PE) ;
  
    print_mat( Heff, " H effective ") ;
    h_vec.push_back( Heff) ;

    *(HCI.data() + i) = PE ;
    *(SCI.data() + i) = PO ;
    }

  print_mat( HCI, " CI Hamiltonian ") ;
  print_mat( SCI, " CI Overlap ") ;

  H_KCI.compute( HCI, SCI) ;
  std::cout << " H_KCI.eigenvalues() " << std::endl ;
  std::cout << H_KCI.eigenvalues() << std::endl ;
  std::cout << " H_KCI.eigenvectors() " << std::endl ;
  std::cout << H_KCI.eigenvectors() << std::endl ;
  e_vec = H_KCI.eigenvectors().col( 0) ;
  print_mat( e_vec, " Lowest Eigenvector ") ;
  CI_S = e_vec.adjoint()*SCI*e_vec ;
  std::cout << " CI_S " << CI_S << std::endl ;
  /* Generate the Effective Hamiltonian for the K and N projected state */
  Heff = (std::conj(e_vec(0))*e_vec(0)*h_vec[0] + std::conj(e_vec(1))*e_vec(0)*h_vec[1] + std::conj(e_vec(0))*e_vec(1)*h_vec[2] + std::conj(e_vec(1))*e_vec(1)*h_vec[3])/CI_S ;
  print_mat( Heff, " CI Hamiltonian ") ;
  ring_shiekh_K_projection_time.end() ;

  return ;

} ;

template void ring_shiekh_K( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd, trapezoid*&, double&, diis_control&, int&) ;

template < class matrix>
void Xring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  This is a working version of the routine however inefficient.  DO NOT ALTER
*/
  int iter = 1, iter_N = 0 ;
  int in, inmb = ngrid->ns() ;
  double b_ll, b_ul ;
  double pnum = nele/d2 ;
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

  int iter = 1, dbas = 2*nbas ;
  int in, inmb = ngrid->ns() ;
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

template < class matrix>
void general_derivative_testing( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::VectorXd> eigval, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) {
/*
  Implementing the more general version of the projected HFB effective Hamiltonian.  MEMORY HUNGRY
*/
  int in, inmb = ngrid->ns(), iter = 0, dbas = 2*nbas ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag ;
//  diis<cd> fdiis( 2*nbas, diis_cntrl.ndiisv, diis_cntrl.diistype) ;
  cd nrm, c_junk ;
  cd energy, pphase, nphase ;
  cd s_theta, intE_g, intx_g, fnx ;
  cd etr, gtr, dtr, d_phi ;
  cd shift ;
//
  cd c_p5p5 = cd( d1/d2, d1/d2) ;
//
  std::vector<cd> olap( inmb), x_phi( inmb), w_phi( inmb), e_hf( inmb), e_pr( inmb) ;
  std::vector<Eigen::MatrixXcd> Y_phi( inmb), Y_kap( inmb) ;
  std::vector<Eigen::MatrixXcd> R_phi( inmb), r_phi( inmb), k_phi( inmb), kbar_phi( inmb) ;
  std::vector<Eigen::MatrixXcd> C_phi( inmb), C_phi_i( inmb), f_phi( inmb), G_phi( inmb), D_phi( inmb), Dbar_phi( inmb) ;
  Eigen::MatrixXcd p_rho, p_kappa, t, mu ;
  Eigen::MatrixXcd epsilonN, gammaN, lambdaN, FockN, DeltaN ;
  Eigen::MatrixXcd I, rho_i, kappa_i ;
  Eigen::MatrixXcd int_Y_phi, int_Y_kap ;
  Eigen::MatrixXcd scr1, scr2, scr3, R ;
  Eigen::MatrixXcd F_rho, F_kap, D_rho, D_kap, Heff, H ;

  time_dbg general_derivative_testing_time = time_dbg("General derivative testing") ;

/*
  Eigen Matrices to accumulate quantities
*/

  I.resize( 2*nbas, 2*nbas) ;
  scr1.resize( 2*nbas, 2*nbas) ;
  scr2.resize( 2*nbas, 2*nbas) ;
  scr3.resize( 4*nbas, 4*nbas) ;
  rho_i.resize( 2*nbas, 2*nbas) ;
  kappa_i.resize( 2*nbas, 2*nbas) ;
  int_Y_phi.resize( 2*nbas, 2*nbas) ;
  int_Y_kap.resize( 2*nbas, 2*nbas) ;
  t.resize( 2*nbas, 2*nbas) ;
  p_rho.resize( 2*nbas, 2*nbas) ;
  p_kappa.resize( 2*nbas, 2*nbas) ;
  F_rho.resize( 2*nbas, 2*nbas) ;
  D_rho.resize( 2*nbas, 2*nbas) ;
  F_kap.resize( 2*nbas, 2*nbas) ;
  D_kap.resize( 2*nbas, 2*nbas) ;
  H.resize( 4*nbas, 4*nbas) ;
  Heff.resize( 4*nbas, 4*nbas) ;
  R.resize( 4*nbas, 4*nbas) ;
  mu.resize( 4*nbas, 4*nbas) ;

  I.setIdentity() ;
  mu.setIdentity() ;
  mu.block( 0, 0, 2*nbas, 2*nbas) *= -z1 ;

  rho_i = rho.inverse() ;
  kappa_i = kappa.inverse() ;

  R.block( 0, 0, 2*nbas, 2*nbas) = rho ;
  R.block( 0, 2*nbas, 2*nbas, 2*nbas) = kappa ;
  R.block( 2*nbas, 0, 2*nbas, 2*nbas) = kappa.adjoint() ;
  R.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = I - rho.conjugate() ;

  scr3.block( 0, 0, 2*nbas, 2*nbas) = -rho*kappa_i.conjugate() ;
  scr3.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = rho.conjugate()*kappa_i ;
  scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
  scr3.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;

  /* norm for pfaffians */
  nrm = z1/pfaffian_H( scr3) ;
  std::cout << " pfaffian norm " << nrm << std::endl ;

  do {

     p_rho = rho ;
     p_kappa = kappa ;

    int_Y_phi.setZero() ;
    int_Y_kap.setZero() ;
    F_rho.setZero() ;
    D_rho.setZero() ;
    F_kap.setZero() ;
    D_kap.setZero() ;
    intE_g = z0 ;
    intx_g = z0 ;

    for ( in = 0; in < inmb; in++){
      fnx = static_cast<cd>( ngrid->q()) ;
      std::cout << fnx << std::endl ;
      /* w_phi */
      c_junk = ngrid->w() ;
      w_phi[in] = c_junk ;
      d_phi = std::exp( -zi*fnx*static_cast<cd>(nele))/( z2*zpi) ;
      /* R_phi */
      scr1 = I*std::exp( zi*fnx) ;
      R_phi[in] = scr1 ;
      /* C_phi_i */
      scr1 = -kappa*R_phi[in].conjugate()*kappa.conjugate()*R_phi[in].adjoint() + rho*R_phi[in]*rho*R_phi[in].adjoint() ;
      C_phi_i[in] = scr1 ;
      /* C_phi */
      scr1 = C_phi_i[in].inverse() ;
      C_phi[in] = scr1 ;
      /* r_phi */
      scr1 = R_phi[in]*rho*R_phi[in].adjoint()*C_phi[in]*rho ;
      r_phi[in] = scr1 ;
      /* k_phi */
      scr1 = R_phi[in]*rho*R_phi[in].adjoint()*C_phi[in]*kappa ;
      k_phi[in] = scr1 ;
      /* kbar_phi */
      scr1 = R_phi[in].conjugate()*kappa.conjugate()*R_phi[in].adjoint()*C_phi[in]*rho ;
      kbar_phi[in] = scr1 ;

      scr1 = R_phi[in]*kappa ;
      scr2 = scr1.inverse() ;
      scr3.block( 0, 0, 2*nbas, 2*nbas) = -R_phi[in]*rho*scr2.conjugate() ;
      scr3.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = rho.conjugate()*kappa_i ;
      scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
      scr3.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;

      /* Olap */
      c_junk = nrm*pfaffian_H( scr3) ;

      olap[in] = c_junk ;

      /* x_phi */
      c_junk = d_phi*olap[in] ;
      x_phi[in] = c_junk ;

      /* Y_phi and int dphi Y_phi */
      scr1 = R_phi[in]*rho*R_phi[in].adjoint()*C_phi[in] + R_phi[in].adjoint()*C_phi[in]*rho*R_phi[in] ;
      Y_phi[in] = scr1/z2 ;
      int_Y_phi += w_phi[in]*x_phi[in]*Y_phi[in] ;

      /* Y_kap and int dphi Y_kap */
      scr1 = R_phi[in].adjoint()*C_phi[in]*kappa*R_phi[in].conjugate() ;
      Y_kap[in] = scr1/z2 ;
      int_Y_kap += w_phi[in]*x_phi[in]*Y_kap[in] ;

      /* G_phi & D_phi */
      scr1 = r_phi[in] ;
      scr2 = k_phi[in] ;
      W->contract( scr1, scr2) ;

      scr3 = W->getG() ;
      scr1 = scr3.block( 0, 0, 2*nbas, 2*nbas) ;
      G_phi[in] = scr1 ;
      scr2 = scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
      D_phi[in] = scr2 ;

      scr2 = kbar_phi[in] ;
/* Overkill for now */
      W->contract( scr1, scr2) ;
      scr3 = W->getG() ;
      scr2 = scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) ;

      Dbar_phi[in] = scr2 ;

      scr1 = h + G_phi[in] ;
      f_phi[in] = scr1 ;

      scr2 = h + f_phi[in] ;
      scr1 = scr2*r_phi[in] ;

      c_junk = scr1.trace()/z2 ;
      e_hf[in] = c_junk ;

      scr1 = D_phi[in]*kbar_phi[in].adjoint() + D_phi[in].transpose()*kbar_phi[in].conjugate() ;
      c_junk = -scr1.trace()/z2 ;
      e_pr[in] = c_junk ;

      intx_g += w_phi[in]*x_phi[in] ;
      intE_g += w_phi[in]*x_phi[in]*( e_hf[in] - e_pr[in]/z2) ;
      }

    ngrid->set_s() ;

    for ( in = 0; in < inmb; in++){
      std::cout << e_hf[in] << std::endl ;
      }
    std::cout << " e_pr " << std::endl << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << e_pr[in]/z2 << std::endl ;
      }
    std::cout << std::endl ;
    std::cout << " Weights " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      std::cout << w_phi[in] << std::endl ;
      }
    std::cout << std::endl ;
    std::cout << " x(phi) " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      std::cout << x_phi[in] << std::endl ;
      }
    std::cout << std::endl ;
    std::cout << " Olap " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      std::cout << olap[in] << std::endl ;
      }
    std::cout << std::endl ;
    std::cout << " R_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( R_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " r_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( r_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " k_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( k_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " kbar_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( kbar_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " C_phi_i " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( C_phi_i[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " C_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( C_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " G_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( G_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " D_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( D_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " Dbar_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( Dbar_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " Y_phi " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( Y_phi[in]) ;
      }
    std::cout << std::endl ;
    std::cout << " Y_kappa " << std::endl ;
    for ( in = 0; in < inmb; in++){
      std::cout << " " << in << std::endl ;
      print_mat( Y_kap[in]) ;
      }

    int_Y_phi /= intx_g ;
    int_Y_kap /= intx_g ;

    print_mat( int_Y_phi, "int_Y_phi") ;
    print_mat( int_Y_kap, "int_Y_kap") ;

    for ( in = 0; in < inmb; in++){
/*
      F_rho += w_phi[in]*x_phi[in]*((Y_phi[in] - int_Y_phi)*e_hf[in] + R_phi[in].adjoint()*C_phi[in]*rho*f_phi[in]*( I - r_phi[in])*R_phi[in] + ( I - r_phi[in])*f_phi[in]*R_phi[in]*rho*R_phi[in].adjoint()*C_phi[in])/intx_g ;
      D_rho += w_phi[in]*x_phi[in]*( -(Y_phi[in] - int_Y_phi)*e_pr[in] + (-R_phi[in].adjoint()*C_phi[in]*rho*D_phi[in]*kbar_phi[in].adjoint()*R_phi[in] + r_phi[in]*D_phi[in]*R_phi[in].conjugate()*kappa.conjugate()*R_phi[in].adjoint()*C_phi[in] - D_phi[in]*R_phi[in].conjugate()*kappa.conjugate()*R_phi[in].adjoint()*C_phi[in] - R_phi[in].adjoint()*C_phi[in]*kappa*Dbar_phi[in].adjoint()*r_phi[in]*R_phi[in] + R_phi[in].adjoint()*C_phi[in]*kappa*Dbar_phi[in].adjoint()*R_phi[in] - k_phi[in]*Dbar_phi[in].adjoint()*R_phi[in]*rho*C_phi[in]) )/(intx_g*z2) ;
      F_kap += w_phi[in]*x_phi[in]*( ( -Y_kap[in]/z2 + int_Y_kap/z2)*e_hf[in] + R_phi[in].adjoint()*C_phi[in]*rho*f_phi[in]*k_phi[in]*R_phi[in].conjugate())/intx_g ;
      D_kap += w_phi[in]*x_phi[in]*( ( -Y_kap[in]/z2 + int_Y_kap/z2)*e_pr[in] + R_phi[in].adjoint()*C_phi[in]*rho*D_phi[in]*r_phi[in].transpose()*R_phi[in].conjugate() - R_phi[in].adjoint()*C_phi[in]*kappa*Dbar_phi[in].adjoint()*k_phi[in]*R_phi[in].conjugate())/intx_g ;
*/
      D_rho = -R_phi[in].adjoint()*C_phi[in]*rho*D_phi[in]*kbar_phi[in].adjoint()*R_phi[in] ;
      print_mat( D_rho, "R C p D kb R") ;
      D_rho = (I - r_phi[in])*D_phi[in]*R_phi[in].conjugate()*kappa.conjugate()*R_phi[in].adjoint()*C_phi[in] ;
      print_mat( D_rho, "( 1 - p) D R^* k^* R^t C") ;
      D_rho = R_phi[in].adjoint()*C_phi[in]*kappa*Dbar_phi[in].adjoint()*(I - r_phi[in])*R_phi[in] ;
      print_mat( D_rho, "R^t C k Db ( 1 - p) R") ;
      D_rho = k_phi[in]*Dbar_phi[in].adjoint()*R_phi[in]*rho*R_phi[in].adjoint()*C_phi[in] ;
      print_mat( D_rho, " k Db R p R^t C") ;
      }

    return ;

    Heff.block( 0, 0, 2*nbas, 2*nbas) = (F_rho + F_rho.adjoint())/z2 + (D_rho + D_rho.adjoint())/z2 ;
    Heff.block( 0, 2*nbas, 2*nbas, 2*nbas) = (-F_kap + F_kap.transpose())/z2 + (D_kap - D_kap.transpose())/z2  ;
    Heff.block( 2*nbas, 0, 2*nbas, 2*nbas) = -Heff.block( 0, 2*nbas, 2*nbas, 2*nbas).conjugate() ;
    Heff.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -Heff.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;

    /* Some crude level shifting */
    scr3.setIdentity() ;
    Heff += -z4*R ;
    Heff += z4*( scr3 - R) ;

    H_diag.compute( Heff) ;
    scr3 = H_diag.eigenvectors() ;
    R = scr3.block( 0, 0, 4*nbas, 2*nbas)*scr3.block( 0, 0, 4*nbas, 2*nbas).adjoint() ;

//    chemical_potential( nele, dbas, lambda, H_diag, R, mu, scr3, Heff, rho, H) ;

    rho = (R.block( 0, 0, 2*nbas, 2*nbas) + R.block( 0, 0, 2*nbas, 2*nbas).adjoint())/z2 ;
    kappa = (R.block( 0, 2*nbas, 2*nbas, 2*nbas) - R.block( 0, 2*nbas, 2*nbas, 2*nbas).transpose())/z2 ;

    t = rho - p_rho ;
    energy = t.norm() ;
    t = kappa - p_kappa ;
    energy += t.norm() ;
    std::cout << " Rmsd in denstiy error " << energy << std::endl ;

   if ( energy.real() < 1.0e-7) {
     std::cout << " Converged in the density after " << iter << " iterations. "  << std::endl ;
     std::cout << " Chemical potential: " << lambda << std::endl ;
     std::cout << " Energy " << intE_g << std::endl ;
     break ;
     }

  } while ( ++iter < 300) ;

  general_derivative_testing_time.end() ;

  return ;

} ;

template void general_derivative_testing( int&, Eigen::Ref<Eigen::MatrixXcd>, nbodyint<Eigen::MatrixXcd>*, Eigen::Ref<Eigen::VectorXd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, double&, cd, trapezoid*&, double&, diis_control&, int&) ;

template < class matrix>
void projected_wavefunction_energy( const int& nele, const int& nbas, nbodyint<matrix>* W, trapezoid*& ngrid, const Eigen::Ref<Eigen::MatrixXcd> h, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, const Eigen::Ref<Eigen::MatrixXcd> kappa_i, Eigen::Ref<Eigen::MatrixXcd> R_phi, matrix& r_phi, matrix& k_phi, Eigen::Ref<Eigen::MatrixXcd> kbar_phi, Eigen::Ref<Eigen::MatrixXcd> I, Eigen::Ref<Eigen::MatrixXcd> scr1, Eigen::Ref<Eigen::MatrixXcd> scr2, Eigen::Ref<Eigen::MatrixXcd> scr3, cd& intx_g, cd& intE_g) {
/*
  This is a wrapper for convienence and readability and so most everything is passed preallocated.
  Given a density, pairing matrix, two electron interaction, and integration grids, calculate the energy and overlap.
*/

  int in, inmb = ngrid->ns() ;
  cd fnx, w_phi, x_phi, d_phi, pphase, nphase, olap, nrm, etr ;

  time_dbg projected_wavefunction_energy_time = time_dbg("projected wavefunction energy") ;

  intx_g = z0 ;
  intE_g = z0 ;
  nrm = z0 ;

  ngrid->set_s() ;

  for ( in = 0; in < inmb; in++){
    fnx = static_cast<cd>( ngrid->q()) ;
    w_phi = ngrid->w() ;
    d_phi = std::exp( -zi*fnx*static_cast<cd>( nele))/( z2*zpi) ;
    pphase = std::exp( z2*zi*fnx) ;
    nphase = std::exp( -z2*zi*fnx) ;

    /*
      Build our rotated quantities
    */

    scr1 = I + rho*( pphase - z1) ;
    scr2 = scr1.inverse() ;
    scr1 = pphase*scr2 ;
    r_phi = scr1*rho ;
    k_phi = scr1*kappa ;
    kbar_phi = pphase*kappa*scr1.conjugate() ;
/* This is super janky right now but I will clean it later */
    R_phi = I*std::exp( zi*fnx) ;
    scr1 = R_phi*kappa ;
    scr2 = scr1.inverse() ;
    scr3.block( 0, 0, 2*nbas, 2*nbas) = -R_phi*rho*scr2.conjugate() ;
    scr3.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = kappa_i*rho ;
    scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) = -I ;
    scr3.block( 2*nbas, 0, 2*nbas, 2*nbas) = I ;

    if ( in == 0 ){
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
    etr += scr2.trace()/z2 ;
    scr1 = scr3.block( 0, 2*nbas, 2*nbas, 2*nbas) ;
    scr2 = scr1*kbar_phi.adjoint() ;
    etr += scr2.trace()/z4 ;
    scr2 = scr1.transpose()*kbar_phi.conjugate() ;
    etr += scr2.trace()/z4 ;

/*
  Accumulate the overlap and projected energy
*/
    intE_g += w_phi*x_phi*etr ;
    }

  ngrid->set_s() ;

  projected_wavefunction_energy_time.end() ;

  return ;

} ;

template void projected_wavefunction_energy( const int&, const int&, nbodyint<Eigen::MatrixXcd>*, trapezoid*&, const Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, const Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, cd&, cd&) ;

void chemical_potential ( const int& nele, const int& nbas, double& lambda, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>& H_diag, Eigen::Ref<Eigen::MatrixXcd> R, Eigen::Ref<Eigen::MatrixXcd> mu, Eigen::Ref<Eigen::MatrixXcd> eigvec, const Eigen::Ref<Eigen::MatrixXcd> Heff, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> H) {
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
