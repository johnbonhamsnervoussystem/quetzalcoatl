#include <cmath>
#include <complex>
#include "common.h"
#include "constants.h"
#include "diis.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "evalm.h"
#include "hfrout.h"
#include "init_job.h"
#include "guess.h"
#include "nbodyint.h"
#include <iostream>
#include "qtzcntrl.h"
#include "qtzio.h"
#include "wfn.h"
#include "r12.h"
#include "solver.h"
#include <string>
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include <vector>

/* 
  This set of routines finds a wavefunction by repeated diagonalization of the fock matrix for 
  all the flavors of HF.
*/
void scf_drv( common& com) {

/* 
  Determine the symmetry of the model
    1 - real restricted
    2 - complex restricted
    3 - real unrestricted
    4 - complex unrestricted
    5 - real generalized
    6 - complex generalized

  Determine the type of wavefunction
    10 - Slater Determinant
    20 - Hartree-Fock-Bogoliubov 
*/
  int wfntyp = ( com.methd() / 10) % 2 ;
  int wfnksy = ( com.methd() % 10) % 2 ;
  int wfnssy = ( ( com.methd() % 10) + 1) / 2 ;
  int wfnguess ;
  
  if ( wfntyp == 0 ) {

/*
  HFB wavefunctions
*/

/*
  Generate an inital guess from a converged HF calculation
*/
    wfnguess = 1 ;
    real_SlaDet( com, wfnguess) ;

    if ( wfnksy == 1 ) {

      /* Do real scf */
      real_HFB( com, wfnssy) ;
    } else if ( wfnksy == 0 ) {

      /* Do Complex scf */
      cplx_HFB( com, wfnssy) ;
      }
  } else if (  wfntyp == 1) {

/*
  Hartree-Fock wavefunctions
*/
    if ( wfnksy == 1 ) {

    /* Do real scf */
      real_SlaDet( com, wfnssy) ;
    } else if ( wfnksy == 0 ) {

    /* Do complex scf */
      cplx_SlaDet( com, wfnssy) ;
      }
    } else {

      std::cout << " method = " << com.methd() << std::endl ;
      qtzcntrl::shutdown( "Unrecognized option in scf_drv") ;
      }

  return ;

  } ;

void real_HFB( common& com, int& opt) {
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int nalp = com.nalp() ;
    int maxit_scf = com.mxscfit() ;
    int maxit_pn = com.mxpnit() ;
//    double lambda = com.mu() ;
    double lambda = d0 ;
    double thresh = com.scfthresh() ;
    double lshift = com.lvlshft() ;
    nbodyint<Eigen::MatrixXd>* X ;
    Eigen::MatrixXd h, xs, xsi, p, k, u, v ;
    std::vector<tei>* r12int ;
    wfn< double, Eigen::Dynamic, Eigen::Dynamic> w ;
    time_dbg real_HFB_time = time_dbg("real_HFB") ;

    if ( opt == 1) {

/*
  Restricted HFB
*/
      initialize( 2, opt, com.hamil(), com, h, X, r12int, nbas) ;
      /* 
        Unfortunately generating the initial guess requires an additional call to getXS...
        Something that will have to be cleaned later.  But only for molecular hamiltonians.
        Siggghhh. 
      */
      p.resize( nbas, nbas) ;
      k.resize( nbas, nbas) ;
      p.setZero() ;
      k.setZero() ;
//      xs.resize( nbas, nbas) ;
//      xs = com.getXS() ;
      thermal_guess( nalp, nbas, p, k) ;
//      xsi.resize( nbas, nbas) ;
//      xsi = xs.inverse() ;
//      transform( 2, xsi, p) ;
//      transform( 2, xsi, k) ;
//      xsi.resize( 0, 0) ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      w.eig.resize( 2*nbas) ;

      w.e_scf = rhfbdia( h, X, nbas, nele, p, k, w.moc, w.eig, lambda, lshift, maxit_scf, maxit_pn, thresh) ;

      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, "MO Eigenvalues : ") ;
      p = w.moc.block( nbas, nbas, nbas, nbas) ;
      k = w.moc.block( 0, nbas, nbas, nbas) ;
      print_mat( p, " V ") ;
      print_mat( k, " U ") ;
    } else if ( opt == 2) {

/*
  Unrestricted HFB
*/
      qtzcntrl::shutdown("Unrestricted HFB is not yet implemented") ;
    } else if ( opt == 3) {

/*
  Generalized HFB
*/
      initialize( 2, opt, com.hamil(), com, h, X, r12int, nbas) ;

/*
  Build rho and kappa from a previous HF calculation
*/
      p.resize( 2*nbas, 2*nbas) ;
      k.resize( 2*nbas, 2*nbas) ;
      p.setZero() ;
      k.setZero() ;
      thermal_guess( nalp, nbas, p, k) ;

      w.moc.resize( 4*nbas, 4*nbas) ;
      w.eig.resize( 4*nbas) ;

      w.e_scf = ghfbdia( h, X, nbas, nele, p, k, w.moc, w.eig, lambda, lshift, maxit_scf, maxit_pn, thresh) ;

      com.mu( std::real(lambda)) ;
      std::cout << "Mean Field Energy + NN : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
    } else {
      qtzcntrl::shutdown("Unrecognized option in real_HFB") ;
      }

//    w.wfntyp = opt ;
    w.eig.resize( 0) ;
    w.moc.resize( 0, 0) ;
    xs.resize( 0, 0) ;
    k.resize( 0, 0) ;
    p.resize( 0, 0) ;
    h.resize( 0, 0) ;
    real_HFB_time.end() ;

    return ;

} ;

void real_SlaDet( common& com, int& opt){
/*
    Do scf for real wavefunctions.
      The dimension is done here but inside the routine it takes
      care of other speficic details of the model.
    rhfdia - restricted hartree fock -> 1
    uhfdia - unrestricted hartree fock -> 2
    ghfdia - generalized hartree fock -> 3
*/
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int nalp = com.nalp() ;
    int nbet = com.nbet() ;
    int maxit = com.mxscfit() ;
    double thresh = com.scfthresh() ;
    double lvlshft = com.lvlshft() ;
    bool lshift = com.use_shift() ;
    nbodyint<Eigen::MatrixXd>* X ;
    std::vector<tei>* r12int ;
    Eigen::MatrixXd h, xs, ca, cb ;
    wfn< double, Eigen::Dynamic, Eigen::Dynamic> w ;
    time_dbg real_scf_time = time_dbg("real_scf") ;

    if ( opt == 1 ) {

      initialize( 1, opt, com.hamil(), com, h, X, r12int, nbas) ;
      w.moc.resize( nbas, nbas) ;
      w.moc.setRandom() ;
      w.moc *= d1/d10 ;

      w.eig.resize( nbas) ;

      w.e_scf = rhfdia( h, X, nbas, nele, w.moc, w.eig, lshift, lvlshft, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf << std::endl ;
      std::cout << "Nuclear Repulsion : " << com.nrep() << std::endl ;
      std::cout << "Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, "MO Eigenvalues : ") ;
      print_mat( w.moc, "MO coefficients : ") ;

    } else if ( opt == 2) {

      h.resize( nbas, nbas) ;
      xs.resize( nbas, nbas) ;
      h = com.getH() ;
      xs = com.getXS() ;
      transform( 2, xs, h) ;
      nbodyint<Eigen::MatrixXd>* X = new r12<Eigen::MatrixXd>( r12int, xs, 2, nbas) ;
      ca.resize( nbas, nbas) ;
      cb.resize( nbas, nbas) ;
      ca.setRandom() ;
      cb.setRandom() ;
      ca *= d1/d10 ;
      cb *= d1/d10 ;
      w.eig.resize( 2*nbas) ;

      w.e_scf = uhfdia( h, X, nbas, nalp, nbet, ca, cb, w.eig, maxit, thresh) ;

      w.moc.resize( nbas, 2*nbas) ;
      w.moc.block(0, 0, nbas, nbas) = ca ;
      w.moc.block(0, nbas, nbas, nbas) = cb ;
      cb.resize( 0, 0) ;
      ca.resize( 0, 0) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, " MO Eigenvalues (a/b) : " ) ;
      print_mat( w.moc, " MO Coeffcients (a/b): " ) ;

    } else if ( opt == 3 ) {

      initialize( 1, opt, com.hamil(), com, h, X, r12int, nbas) ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      w.moc.setRandom() ;
      w.moc *= d1/d10 ;

      w.eig.resize( 2*nbas) ;

      w.e_scf = ghfdia( h, X, nbas, nele, w.moc, w.eig, lshift, lvlshft, maxit, thresh) ;

      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, "MO Eigenvalues : ") ;
      print_mat( w.moc, "MO coefficients : ") ;

    } else {
      std::cout << " opt = " << opt << std::endl ;
      qtzcntrl::shutdown("Unrecognized option in real_SlaDet") ;
      }

    w.wfntyp = opt ;
    save_wfn(w) ;
    w.moc.resize( 0, 0) ;
    w.eig.resize( 0) ;
    xs.resize( 0, 0) ;
    h.resize( 0, 0) ;

    real_scf_time.end() ;

    return ;

  } ;

void cplx_HFB( common& com, int& opt) {
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int nalp = com.nalp() ;
    int maxit_scf = com.mxscfit() ;
    int maxit_pn = com.mxpnit() ;
    double thresh = com.scfthresh() ;
    cd lambda = cd ( com.mu(), d0) ;
    cd lshift = -z4 ;
    Eigen::MatrixXcd h, xs, xsi, p, k ;
    std::vector<tei>* r12int ;
    nbodyint<Eigen::MatrixXcd>* X ;
    com.getr12( r12int) ;
    wfn< cd, Eigen::Dynamic, Eigen::Dynamic> w ;
    time_dbg cplx_HFB_time = time_dbg("cplx_HFB") ;

    if ( opt == 1 ) {

/*
  Restricted HFB
*/
      h.resize( 2*nbas, 2*nbas) ;
      p.resize( nbas, nbas) ;
      k.resize( nbas, nbas) ;
      xs.resize( nbas, nbas) ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      w.eig.resize( 2*nbas) ;
      h.setZero() ;
      p.setZero() ;
      k.setZero() ;
      xs.setZero() ;
      xs.real() = com.getXS() ;
      h.block( 0, 0, nbas, nbas).real() = com.getH() ;
      transform( 0, xs, h.block( 0, 0, nbas, nbas)) ;
      h.block( nbas, nbas, nbas, nbas) = -h.block( 0, 0, nbas, nbas) ;
      nbodyint<Eigen::MatrixXcd>* X = new r12<Eigen::MatrixXcd>( r12int, xs, 4, nbas) ;
      thermal_guess( nalp, nbas, p, k) ;
      xsi.resize( nbas, nbas) ;
      xsi = xs.inverse() ;
      transform( 2, xsi, p) ;
      transform( 2, xsi, k) ;
      xsi.resize( 0, 0) ;

      w.e_scf = rhfbdia( h, X, nbas, nele, p, k, w.moc, w.eig, lambda, lshift, maxit_scf, maxit_pn, thresh) ;

      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, "MO Eigenvalues : ") ;
      print_mat( w.moc, "MO coefficients : ") ;

    } else if ( opt == 2 ) {
      qtzcntrl::shutdown("Not yet implemented") ;
    } else if ( opt == 3 ) {

      initialize( 2, opt, com.hamil(), com, h, X, r12int, nbas) ;
      h.resize( 4*nbas, 4*nbas) ;
      xs.resize( 2*nbas, 2*nbas) ;
      p.resize( 2*nbas, 2*nbas) ;
      k.resize( 2*nbas, 2*nbas) ;
      w.moc.resize( 4*nbas, 4*nbas) ;
      w.eig.resize( 4*nbas) ;
      h.setZero() ;
      p.setZero() ;
      k.setZero() ;
      xs.setZero() ;
      xs.block( 0, 0, nbas, nbas).real() = com.getXS() ;
      xs.block( nbas, nbas, nbas, nbas) = xs.block( 0, 0, nbas, nbas) ;
      h.block( 0, 0, nbas, nbas).real() = com.getH() ;
      h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
      transform( 2, xs, h.block( 0, 0, 2*nbas, 2*nbas)) ;
      h.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -h.block( 0, 0, 2*nbas, 2*nbas) ;
      nbodyint<Eigen::MatrixXcd>* X = new r12<Eigen::MatrixXcd>( r12int, xs, 6, nbas) ;

/*
  Build rho and kappa from a previous HF calculation
*/

      thermal_guess( nalp, nbas, p, k) ;
      xsi.resize( 2*nbas, 2*nbas) ;
      xsi = xs.inverse() ;
      transform( 2, xsi, p) ;
      transform( 2, xsi, k) ;
      xsi.resize( 0, 0) ;

      w.e_scf = ghfbdia( h, X, nbas, nele, p, k, w.moc, w.eig, lambda, lshift, maxit_scf, maxit_pn, thresh) ;

      com.mu( std::real(lambda)) ;
      std::cout << "Mean Field Energy + NN : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
    } else {
      qtzcntrl::shutdown("Unrecognized option in cplx_HFB") ;
      }

//    w.wfntyp = wt/2 ;
    save_wfn( w) ;
    w.eig.resize( 0) ;
    w.moc.resize( 0, 0) ;
    k.resize( 0, 0) ;
    p.resize( 0, 0) ;
    xs.resize( 0, 0) ;
    h.resize( 0, 0) ;

    cplx_HFB_time.end() ;

    return ;

} ;

void cplx_SlaDet( common& com, int& opt){
/*  
    Do scf for real wavefunctions.  
      The dimension is done here but inside the routine it takes 
      care of other speficic details of the model.
    rhfdia - restricted hartree fock -> 1
    uhfdia - unrestricted hartree fock -> 2
    ghfdia - generalized hartree fock -> 3 
*/
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int nalp = com.nalp() ;
    int nbet = com.nbet() ;
    int maxit = com.mxscfit() ;
    double thresh = com.scfthresh() ;
    double lvlshft = com.lvlshft() ;
    bool lshift = com.use_shift() ;
    nbodyint<Eigen::MatrixXcd>* X ;
    std::vector<tei>* r12int ;
    Eigen::MatrixXcd h, xs, ca, cb, hfc ;
    wfn< cd, Eigen::Dynamic, Eigen::Dynamic> w ;
    time_dbg cplx_hf_time = time_dbg("complex_hartree_fock") ;

    if ( opt == 1 ) {

      initialize( 1, opt, com.hamil(), com, h, X, r12int, nbas) ;
      w.moc.resize( nbas, nbas) ;
      w.moc.setRandom() ;
      w.moc *= z1/z10 ;

      w.eig.resize( nbas) ;

      w.e_scf = rhfdia( h, X, nbas, nele, w.moc, w.eig, lshift, lvlshft, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf << std::endl ;
      std::cout << "Nuclear Repulsion : " << com.nrep() << std::endl ;
      std::cout << "Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, "MO Eigenvalues : ") ;
      print_mat( w.moc, "MO coefficients : ") ;

    } else if ( opt == 2 ) {

      initialize( 1, opt, com.hamil(), com, h, X, r12int, nbas) ;
      ca.resize( nbas, nbas) ;
      cb.resize( nbas, nbas) ;
      ca.setRandom() ;
      cb.setRandom() ;
      ca *= z1/z10 ;
      cb *= z1/z10 ;
      w.eig.resize( 2*nbas) ;
      w.moc.setZero() ;

      w.e_scf = uhfdia( h, X, nbas, nalp, nbet, ca, cb, w.eig, maxit, thresh) ;

      w.moc.resize( nbas, 2*nbas) ;
      w.moc.block(0, 0, nbas, nbas) = ca ;
      w.moc.block(0, nbas, nbas, nbas) = cb ;
      cb.resize( 0, 0) ;
      ca.resize( 0, 0) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, " MO Eigenvalues ( a/b) : " ) ;
      print_mat( w.moc, " MO Coeffcients ( a/b): " ) ;

    } else if ( opt == 3 ) {

      initialize( 1, opt, com.hamil(), com, h, X, r12int, nbas) ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      w.moc.setRandom() ;

      w.eig.resize( 2*nbas) ;

/*
  Generate a fermi contact term to try and help break symmetry
    - For tomorrow
        Add a poiter to the basis set in common.

*/

//      w.e_scf = ghfdia( h, X, nbas, nele, w.moc, w.eig, lshift, lvlshft, maxit, thresh) ;

      hfc = com.getFC() ;
      w.e_scf = ghfdia_fc( h, hfc, X, nbas, nele, w.moc, w.eig, lshift, lvlshft, maxit, thresh) ;

      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      print_mat( w.eig, "MO Eigenvalues : ") ;
      print_mat( w.moc, "MO coefficients : ") ;

    } else {
      std::cout << " opt = " << opt << std::endl ;
      qtzcntrl::shutdown("Unrecognized option in cplx_SlaDet") ;
      }

    w.wfntyp = opt ;
    save_wfn(w) ;
    w.eig.resize( 0) ;
    w.moc.resize( 0, 0) ;
    xs.resize( 0, 0) ;
    h.resize( 0, 0) ;

    cplx_hf_time.end() ;

    return ;

  } ;

template < class matrix>
double rhfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nele, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, bool lshift, double lvlshft, const int& maxit, const double& thresh) {

  /* 
    Restricted Hartree-Fock solved by repeated diagonalization.
  */
  matrix f, g, p, p_prev ;
  Eigen::SelfAdjointEigenSolver<matrix> f_diag ;
  int iter=0 ;
  int occ, diis_on = 0 ;
  typename matrix::Scalar energy, prev_energy, t, tx, shift, zero, two ;
  diis<typename matrix::Scalar> p_diis ( nbasis*nbasis, 6) ;
  time_dbg rhfdia_time = time_dbg("rhfdia") ;

  zero = static_cast<typename matrix::Scalar>( d0) ;
  two = static_cast<typename matrix::Scalar>( d2) ;
  energy = zero ;
  occ = nele/2 ;
  if ( lshift){
    shift = static_cast<typename matrix::Scalar>( lvlshft) ;
    }
  f.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p.resize( nbasis, nbasis) ;
  p_prev.resize( nbasis, nbasis) ;

  /* At some point its probably better just to pass in the density */
  p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;

  while ( iter++ < maxit ) {
    prev_energy = energy ;
    p_prev = p ;
    W->contract( p) ;
    g = W->getG() ;
    f = h + g ;
    g = p*( h + f) ;

    energy = g.trace() ;

    if ( false ){
      f += -shift*p ;
      }
    f_diag.compute( f) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;

/*

    if ( diis_on == 0){
      tx = energy - prev_energy ;
      t = std::abs( tx) ;
      std::cout << " Energy difference between iterations is :  " << t << std::endl ;
      if ( std::real(t) < 1.0e-3){
        std::cout << " Turning on DIIS " << std::endl ;
        diis_on = 1 ;
        }
      }

    if ( diis_on == 1){
      p_diis.update( p, diis_on) ;
      }
*/

    g = p - p_prev ;
    t = g.norm() ;
    std::cout << "  rms difference in the densities: " << t << std::endl ;
    if ( std::real(t) < thresh ) { break ;}

    }

  std::cout << " Number of iterations : " << iter << std::endl ;
  p = two*c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
  print_mat( p, " Final Density") ;

  eig = f_diag.eigenvalues() ;
  p_prev.resize( 0, 0) ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  rhfdia_time.end() ;

  return std::real(energy) ;

} ;

template double rhfdia( const Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*, const int&, const int&, Eigen::MatrixXd&, Eigen::Ref<Eigen::VectorXd>, bool, double, const int&, const double&) ;

template double rhfdia( const Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*, const int&, const int&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::VectorXd>, bool, double, const int&, const double&) ;

template < class matrix>
double uhfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nalp, const int& nbet, matrix& c_a, matrix& c_b, Eigen::Ref<Eigen::VectorXd> eig, const int& maxit, const double& thresh){

  /* Real unrestricted Hartree-Fock solved by repeated diagonalization. */
  matrix f_a, f_b, p_t, p_a, p_b, g, p_pa, p_pb ;
  Eigen::SelfAdjointEigenSolver<matrix> f_diag ;

  int iter=0 ;
  typename matrix::Scalar energy, t ;
  time_dbg ruhfdia_time = time_dbg("uhfdia") ;

  f_a.resize( nbasis, nbasis) ;
  f_b.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p_t.resize( nbasis, nbasis) ;
  p_a.resize( nbasis, nbasis) ;
  p_b.resize( nbasis, nbasis) ;
  p_pa.resize( nbasis, nbasis) ;
  p_pb.resize( nbasis, nbasis) ;

  /* If c_a and c_b have values use them as the initial guess */
  p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
  p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;

  while ( iter++ < maxit ) {
    p_pa = p_a ;
    p_pb = p_b ;
    p_t = p_a + p_b ;
    W->contract( p_t, p_a) ;
    g = W->getG() ;
    f_a = h + g ;
    g = p_a*f_a ;
    energy = g.trace() ;
    W->contract( p_t, p_b) ;
    g = W->getG() ;
    f_b = h + g ;
    g = p_b*f_b ;
    energy += g.trace() ;
    g = p_t*h ;
    energy += g.trace() ;

    f_diag.compute( f_a) ;
    c_a = f_diag.eigenvectors() ;
    f_diag.compute( f_b) ;
    c_b = f_diag.eigenvectors() ;
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;

    g = p_a + p_b - p_pa - p_pb ;
    t = g.norm() ;
    if ( std::real(t) < thresh ) { break ;}
    }

  std::cout << " Number of iterations : " << iter << std::endl ;

  eig.head(nbasis) = f_diag.eigenvalues() ;
  eig.tail(nbasis) = f_diag.eigenvalues() ;
  p_t.resize( 0, 0) ;
  p_pb.resize( 0, 0) ;
  p_pa.resize( 0, 0) ;
  p_b.resize( 0, 0) ;
  p_a.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f_b.resize( 0, 0) ;
  f_a.resize( 0, 0) ;

  ruhfdia_time.end() ;

  return std::real(energy)/d2 ;

} ;

template double uhfdia( const Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*, const int&, const int&, const int&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::Ref<Eigen::VectorXd>, const int&, const double&) ;

template double uhfdia( const Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*, const int&, const int&, const int&, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::VectorXd>, const int&, const double&) ;

template < class matrix>
double ghfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nele, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, bool lshift, double lvlshft, const int& maxit, const double& thresh) {

/* 
  Generalized Hartree-Fock solved by repeated diagonalization. 
*/
  matrix f, g, p, p_prev ;
  Eigen::SelfAdjointEigenSolver<matrix> f_diag ;
  int iter=0 ;
  int nbas, diis_on = 0 ;
  typename matrix::Scalar energy, prev_energy, t, tx, shift, zero ;
  diis<typename matrix::Scalar> p_diis ( 2*2*nbasis*nbasis, 6) ;
  time_dbg ghfdia_time = time_dbg("ghfdia") ;

  zero = static_cast<typename matrix::Scalar>( d0) ;
  energy = zero ;
  if ( lshift){
    shift = static_cast<typename matrix::Scalar>(lvlshft) ;
    }
  nbas = nbasis*2 ;
  f.resize( nbas, nbas) ;
  g.resize( nbas, nbas) ;
  p.resize( nbas, nbas) ;
  p_prev.resize( nbas, nbas) ;

  /* Change this at some point */
  p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;

  while ( iter++ < maxit ) {
    prev_energy = energy ;
    p_prev = p ;
    W->contract( p) ;
    g = W->getG() ;
    f = h + g ;
    g = p*( h + f) ;

    energy = g.trace() ;
    if ( false ){
      f += -shift*p ;
      }
    f_diag.compute( f) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;

    if ( diis_on == 0){
      tx = energy - prev_energy ;
      t = std::abs( tx) ;
      std::cout << " Energy difference between iterations is :  " << t << std::endl ;
      if ( std::real(t) < 1.0e-3){
        std::cout << " Turning on DIIS " << std::endl ;
        diis_on = 1 ;
        }
      }
    if ( diis_on == 1){
      p_diis.update( p, diis_on) ;
      }

    g = p - p_prev ;
    t = g.norm() ;
    std::cout << "  rms difference in the densities: " << t << std::endl ;
    if ( std::real(t) < thresh ) { break ;}

    }

  std::cout << " Number of iterations : " << iter << std::endl ;
  p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
  print_mat( p, " Final Density") ;

  eig = f_diag.eigenvalues() ;
  p_prev.resize( 0, 0) ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  ghfdia_time.end() ;

  return std::real(energy)/d2 ;

} ;

template double ghfdia( const Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*, const int&, const int&, Eigen::MatrixXd&, Eigen::Ref<Eigen::VectorXd>, bool, double, const int&, const double&) ;

template double ghfdia( const Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*, const int&, const int&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::VectorXd>, bool, double, const int&, const double&) ;

double ghfdia_fc( const Eigen::Ref<Eigen::MatrixXcd> h, const Eigen::Ref<Eigen::MatrixXcd> fc, nbodyint<Eigen::MatrixXcd>* W, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, bool lshift, double lvlshft, const int& maxit, const double& thresh) {

/* 
  Generalized Hartree-Fock solved by repeated diagonalization. 
*/
  Eigen::MatrixXcd f, g, p, p_prev ;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> f_diag ;
  int iter=1 ;
  int nbas, diis_on = 0 ;
  cd energy, prev_energy, t, tx, shift, zero ;
  diis<cd> p_diis ( 2*2*nbasis*nbasis, 6) ;
  time_dbg ghfdia_fc_time = time_dbg("ghfdia_fc") ;

  energy = z0 ;
  if ( lshift){
    shift = static_cast<cd>(lvlshft) ;
    }
  nbas = nbasis*2 ;
  f.resize( nbas, nbas) ;
  g.resize( nbas, nbas) ;
  p.resize( nbas, nbas) ;
  p_prev.resize( nbas, nbas) ;

  /* Change this at some point */
  p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;

  while ( iter++ < maxit ) {
    prev_energy = energy ;
    p_prev = p ;
    W->contract( p) ;
    g = W->getG() ;
    f = h + g ;
    g = p*( h + f) ;

    energy = g.trace() ;
    if ( iter < 15 && iter % 2 == 0 ){
      f += fc/z4 ;
      }
    f_diag.compute( f) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;

    if ( diis_on == 0){
      tx = energy - prev_energy ;
      t = std::abs( tx) ;
      std::cout << " Energy difference between iterations is :  " << t << std::endl ;
      if ( std::real(t) < 1.0e-3){
        std::cout << " Turning on DIIS " << std::endl ;
        diis_on = 1 ;
        }
      }
    if ( diis_on == 1){
      p_diis.update( p, diis_on) ;
      }

    g = p - p_prev ;
    t = g.norm() ;
    std::cout << "  rms difference in the densities: " << t << std::endl ;
    if ( std::real(t) < thresh ) { break ;}

    }

  std::cout << " Number of iterations : " << iter << std::endl ;
  p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
  print_mat( p, " Final Density") ;

  eig = f_diag.eigenvalues() ;
  p_prev.resize( 0, 0) ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  ghfdia_fc_time.end() ;

  return std::real(energy)/d2 ;

} ;

template < class matrix, class z>
double rhfbdia( const matrix& h, nbodyint<matrix>* X, const int& nbasis, const int& nele, matrix& p, matrix& k, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, z& lambda, z& lvlshift, const int& maxit_scf, const int& maxit_pn, const double& thresh) {
/* 
  Restricted Hartree-Fock Bogoliubov

  - There are two convergence criteria here to achieve self-consistency. 

    For each density, the chemical potential is adjusted until the 
    preferred number of particles is found.  Once the number of particles
    is correct, the density is recomputed and checked against the previous
    iteration.  If it is below the convergence criteria then we are done.

  Input :
    h - core hamiltonian
    X - nbody interaction
    nbasis - the number of spin free basis functions
    nele - the target for particle number in HFB
    p - input guess for the density matrix
    k - input guess for the pairing density matrix
    c - container for the coefficients
    eig - container for the eigenvalues
    lambda - chemical potential.  This take the intial guess
    maxit_scf - max iterations before scf convergence is stopped
    maxit_pn - max iterations before particle number convergence is stopped
    thresh - threshold for convergence of the density

  Local :
    H - HFB Hamiltonian
    p_prev - Stores the previous iteration of the normal density
    k_prev - Stores the previous iteration of the abnormal density
    R - the generalized Density matrix
    G - Self-Consistent Field 
    D - Pairing Field 
    W - Container for two particle interactions
    mu - anti-symmetric overlap for updating chemical potential
    H_diag - Eigensolver object
    t - scratch space
    iter_d - iterations on density convergence
    iter_N - iterations on particle number
    b_ul/b_ll - Bisection method upper and lower limits
    N - number of particles
    lvlshift - level shift parameter

  Return Value :
    energy - final HFB energy/ this also stores the rms error of the density
	in order to check convergence

  To Do :
    Include options to pass in lvlshift
    Include logic to turn lvlshifting on
    Include options to turn on diis
*/

  matrix H ;
  matrix Heff ;
  matrix p_prev ;
  matrix k_prev ;
  matrix R ;
  matrix G ;
  matrix D ;
  matrix W ;
  matrix I ;
  matrix mu ;
  matrix t ;
  Eigen::SelfAdjointEigenSolver<matrix> H_diag ;

  int iter_N = 0, iter_d = 0 ;
  double pnum = static_cast<double>(nele)/d2 ;
  z b_ul, b_ll ;
  typename matrix::Scalar energy, N ;
  time_dbg rhfbdia_time = time_dbg("rhfbdia") ;

/*
  Alloate space for the local matrices
*/

  H.resize( 2*nbasis, 2*nbasis) ;
  Heff.resize( 2*nbasis, 2*nbasis) ;
  p_prev.resize( nbasis, nbasis) ;
  k_prev.resize( nbasis, nbasis) ;
  G.resize( nbasis, nbasis) ;
  D.resize( nbasis, nbasis) ;
  I.resize( nbasis, nbasis) ;
  W.resize( 2*nbasis, 2*nbasis) ;
  R.resize( 2*nbasis, 2*nbasis) ;
  mu.resize( 2*nbasis, 2*nbasis) ;
  t.resize( nbasis, nbasis) ;

  I.setIdentity() ;
  p_prev.setZero() ;
  k_prev.setZero() ;
  R.setZero() ;
  mu.setZero() ;
  mu.block( 0, 0, nbasis, nbasis) = -I ;
  mu.block( nbasis, nbasis, nbasis, nbasis) = I ;

/*
  If no initial guess is given, use the core Hamiltonian
*/

  if ( p.isZero(0)) {
    H = h ;
    H_diag.compute( H) ;
    c = H_diag.eigenvectors() ;
    p = c.block( nbasis, nbasis, nbasis, nbasis).conjugate()*c.block( nbasis, nbasis, nbasis, nbasis).transpose() ;
    k = c.block( nbasis, nbasis, nbasis, nbasis).conjugate()*c.block( 0, nbasis, nbasis, nbasis).transpose() ;
    }

    X->contract( p, k) ;
    W = X->getG() ;
    t = (static_cast<z>(d2)*h.block( 0, 0, nbasis, nbasis) + W.block( 0, 0, nbasis, nbasis))*p ;
    energy = t.trace() ;
    std::cout << " HF Energy " << energy << std::endl ;
    t = k.transpose()*W.block( 0, nbasis, nbasis, nbasis) ;
    std::cout << " Pairing Energy " << t.trace() << std::endl ;
    energy += t.trace() ;
    std::cout << " Total Electronic Energy " << energy << std::endl ;
    R.block( 0, 0, nbasis, nbasis) = p ;
    R.block( nbasis, nbasis, nbasis, nbasis) = I - p.conjugate() ;
    R.block( 0, nbasis, nbasis, nbasis) = k ;
    R.block( nbasis, 0, nbasis, nbasis) = -k.conjugate() ;
    p_prev = p ;
    k_prev = k ;

  while ( iter_d++ < maxit_scf ) {
    std::cout << std::endl << std::endl << "Self-Consistency iteration number: " << iter_d << std::endl << std::endl ;

    Heff = h + W - lvlshift*R ;

    iter_N = 0 ;

/*
   Set and initial guess for the limits of the chemical potential
*/

    b_ul = lambda ;
    b_ll = lambda ;

/*
   Set some initial limits
*/

    do {
      b_ll -= d3 ;
      H = Heff + b_ll*mu ;
      H_diag.compute( H) ;
      c = H_diag.eigenvectors() ;
      p = c.block( nbasis, nbasis, nbasis, nbasis).conjugate()*c.block( nbasis, nbasis, nbasis, nbasis).transpose() ;
      N = p.trace() ;
      } while ( static_cast<double>(std::real(N)) > pnum) ;
  
    std::cout << " Lower limit chemical potential " << b_ll << std::endl ;
  
    do {
      b_ul += d3 ;
      H = Heff + b_ul*mu ;
      H_diag.compute( H) ;
      c = H_diag.eigenvectors() ;
      p = c.block( nbasis, nbasis, nbasis, nbasis).conjugate()*c.block( nbasis, nbasis, nbasis, nbasis).transpose() ;
      N = p.trace() ;
      } while ( static_cast<double>(std::real(N)) < pnum) ;
  
    std::cout << " upper limit chemical potential " << b_ul << std::endl ;
  
    lambda = (b_ul + b_ll)/d2 ;
  
    while ( iter_N++ < maxit_pn) {
      H = Heff + lambda*mu ;
      H_diag.compute( H) ;
      c = H_diag.eigenvectors() ;
      p = c.block( nbasis, nbasis, nbasis, nbasis).conjugate()*c.block( nbasis, nbasis, nbasis, nbasis).transpose() ;
      N = p.trace() ;
      if  ( std::abs(static_cast<double>(std::real(N)) - pnum) < 1.0e-5){
        std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
        std::cout << "    chemical potential: " << lambda << std::endl ;
        std::cout << "    Particle Number : " << d2*std::real(N) << std::endl << std::endl ;
        break ;
      } else {
/*
  Bisection method
*/
        if ( static_cast<double>(std::real(N)) - pnum < d0){
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

    k = c.block( nbasis, nbasis, nbasis, nbasis).conjugate()*c.block( 0, nbasis, nbasis, nbasis).transpose() ;
//    print_mat( k, " kappa") ;
/*
  Now compare densities to check convergence.
*/
//    print_mat( p, " rho after iter : " + std::to_string( iter_d)) ;
//    print_mat( k, " kappa after iter : " + std::to_string( iter_d)) ;
    t = p_prev - p ;
    energy = t.norm() ;
    std::cout << " rms difference in rho: " << energy << std::endl ;
    t = k_prev - k ;
    energy += t.norm() ;
    std::cout << " rms difference in kappa: " << t.norm() << std::endl ;
/*
  Check that the density has converged
*/
    if ( std::real(energy) < thresh ) {
      std::cout << " Converged on the density " << std::endl ;
      break ;
      }

    p_prev = p ;
    k_prev = k ;
    R.block( 0, 0, nbasis, nbasis) = p ;
    R.block( nbasis, nbasis, nbasis, nbasis) = I - p.conjugate() ;
    R.block( 0, nbasis, nbasis, nbasis) = k ;
    R.block( nbasis, 0, nbasis, nbasis) = -k.conjugate() ;
    X->contract( p, k) ;
    W = X->getG() ;
//    print_mat( W, " W @ iter " + std::to_string(iter_d)) ;
  }

/*
  Save the eigenvalues and vectors
*/

  t = (d2*h.block( 0, 0, nbasis, nbasis) + W.block( 0, 0, nbasis, nbasis))*p ;
  energy = t.trace() ;
  std::cout << " Hartree-Fock Energy " << energy << std::endl ;
  t = k.transpose()*W.block( 0, nbasis, nbasis, nbasis) ;
  std::cout << " Pairing Energy " << t.trace() << std::endl ;
  energy += t.trace() ;
  std::cout << " Total Electronic Energy " << energy << std::endl ;
  eig = H_diag.eigenvalues().real() ;
  c = H_diag.eigenvectors() ;

  print_mat( p, "final rho") ;
  print_mat( k, "final kappa") ;

/*
  Clean up the memory
*/

  t.resize( 0, 0) ;
  mu.resize( 0, 0) ;
  R.resize( 0, 0) ;
  W.resize( 0, 0) ;
  I.resize( 0, 0) ;
  D.resize( 0, 0) ;
  G.resize( 0, 0) ;
  k_prev.resize( 0, 0) ;
  p_prev.resize( 0, 0) ;
  H.resize( 0, 0) ;

  rhfbdia_time.end() ;

  return std::real(energy) ;

} ;

template double rhfbdia( const Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*, const int&, const int&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::Ref<Eigen::VectorXd>, double&, double&, const int&, const int&, const double&) ;

template double rhfbdia( const Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*, const int&, const int&, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::VectorXd>, cd&, cd&, const int&, const int&, const double&) ;

template < class matrix, class z>
double ghfbdia( const matrix& h, nbodyint<matrix>* X, const int& nbasis, const int& nele, matrix& p, matrix& k, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, z& lambda, z& lvlshift, const int& maxit_scf, const int& maxit_pn, const double& thresh) {
/* 
  Generalized Hartree-Fock Bogoliubov

  - There are two convergence criteria here to achieve self-consistency. 

    For each density, the chemical potential is adjusted until the 
    preferred number of particles is found.  Once the number of particles
    is correct, the density is recomputed and checked against the previous
    iteration.  If it is below the convergence criteria then we are done.

  Input :
    h - core hamiltonian
    s - overlap
    intarr - two electron integrals
    nbasis - the number of spin free basis functions
    nele - the target for particle number in HFB
    c - container for the coefficients
    eig - container for the eigenvalues
    lambda - chemical potential.  This take the intial guess
    maxit - maxiterations before convergence is stopped
    thresh - threshold for convergence o the density

  Local :
    H - HFB Hamiltonian
    p_prev - Stores the previous iteration of the normal density
    k_prev - Stores the previous iteration of the abnormal density
    R - the top half of the generalized density matrix
    G - Self-Consistent Field 
    D - Pairing Field 
    W - Container for Self-Consistant and Pairing fields
    mu - anti-symmetric overlap for updating chemical potential
    No - R*s for calculating particle number
    H_diag - Eigensolver object
    t - scratch space
    iter_d - iterations on density convergence
    iter_N - iterations on particle number
    b_ul/b_ll - Bisection method upper and lower limits
    N - number of particles
    N_p - number of particles from a previous iteration

  Return Value :
    energy - final HFB energy/ this also stores the rms error of the density
	in order to check convergence

  To DO :
    Add an option to make the pairing field repulsive.
*/

  matrix H ;
  matrix p_prev ;
  matrix k_prev ;
  matrix R ;
  matrix G ;
  matrix D ;
  matrix W ;
  matrix mu ;
  matrix I ;
  matrix t ;
  Eigen::SelfAdjointEigenSolver<matrix> H_diag ;

  int iter_N = 0, iter_d = 0 ;
  double b_ul, b_ll ;
  typename matrix::Scalar energy, N ;
  time_dbg ghfbdia_time = time_dbg("ghfbdia") ;

/*
  Alloate space for the local matrices
*/

  H.resize( 4*nbasis, 4*nbasis) ;
  p_prev.resize( 2*nbasis, 2*nbasis) ;
  k_prev.resize( 2*nbasis, 2*nbasis) ;
  R.resize( 4*nbasis, 4*nbasis) ;
  I.resize( 2*nbasis, 2*nbasis) ;
  G.resize( 2*nbasis, 2*nbasis) ;
  D.resize( 2*nbasis, 2*nbasis) ;
  W.resize( 4*nbasis, 4*nbasis) ;
  mu.resize( 4*nbasis, 4*nbasis) ;
  t.resize( 2*nbasis, 2*nbasis) ;

  p_prev.setZero() ;
  k_prev.setZero() ;
  I.setIdentity() ;
  mu.setZero() ;
  mu.block( 0, 0, 2*nbasis, 2*nbasis) = -I ;
  mu.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis) = I ;

/*
  If no initial guess is given, use the core Hamiltonian
*/

  if ( p.isZero(0)) {
    H = h ;
    H_diag.compute( H) ;
    c = H_diag.eigenvectors() ;
    p = c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).conjugate()*c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).transpose() ;
    k = c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).conjugate()*c.block( 0, 2*nbasis, 2*nbasis, 2*nbasis).transpose() ;
    }

    X->contract( p, k) ;
    W = X->getG() ;
    t = (d2*h.block( 0, 0, 2*nbasis, 2*nbasis) + W.block( 0, 0, 2*nbasis, 2*nbasis))*p ;
    energy = t.trace() ;
    std::cout << " HF Energy " << std::real(energy)/d2 << std::endl ;
    t = k.adjoint()*W.block( 0, 2*nbasis, 2*nbasis, 2*nbasis) ;
    std::cout << " Pairing Energy " << std::real(t.trace())/d2 << std::endl ;
    energy += t.trace() ;
    std::cout << " Total Electronic Energy " << std::real(energy)/d2 << std::endl ;
    R.block( 0, 0, 2*nbasis, 2*nbasis) = p ;
    R.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis) = I - p.conjugate() ;
    R.block( 0, 2*nbasis, 2*nbasis, 2*nbasis) = k ;
    R.block( 2*nbasis, 0, 2*nbasis, 2*nbasis) = -k.conjugate() ;
    p_prev = p ;
    k_prev = k ;

  while ( iter_d++ < maxit_scf ) {
    std::cout << std::endl << std::endl << "Self-Consistency iteration number: " << iter_d << std::endl << std::endl ;
    iter_N = 0 ;
/*
   Set and initial guess for the limits of the chemical potential
*/
    b_ul = std::real(lambda) ;
    b_ll = std::real(lambda) ;

/*
   Set some initial limits
*/

    do {
      b_ll -= d3 ;
      H = h + W + lvlshift*R + static_cast<z>(b_ll)*mu ;
      H_diag.compute( H) ;
      c = H_diag.eigenvectors() ;
      p = c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).conjugate()*c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).transpose() ;
      N = p.trace() ;
      } while ( static_cast<double>(std::real(N)) > static_cast<double>(nele)) ;

    std::cout << " Lower limit chemical potential " << b_ll << std::endl ;

    do {
      b_ul += d3 ;
      H = h + W + lvlshift*R + static_cast<z>(b_ul)*mu ;
      H_diag.compute( H) ;
      c = H_diag.eigenvectors() ;
      p = c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).conjugate()*c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).transpose() ;
      N = p.trace() ;
      } while ( static_cast<double>(std::real(N)) < static_cast<double>(nele)) ;

    std::cout << " upper limit chemical potential " << b_ul << std::endl ;

    lambda = static_cast<z>((b_ul + b_ll)/d2) ;
    H = h + W + lvlshift*R + lambda*mu ;

    while ( iter_N++ < maxit_pn) {
      H_diag.compute( H) ;
      c = H_diag.eigenvectors() ;
      p = c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).conjugate()*c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).transpose() ;
      N = p.trace() ;
      if ( std::abs(static_cast<double>(std::real(N)) - static_cast<double>(nele)) < 1.0e-5){
        std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
        std::cout << "    chemical potential: " << lambda << std::endl ;
        std::cout << "    Particle Number : " << N << std::endl << std::endl ;
        break ;
      } else {

/*
  Bisection method
*/

        if ( static_cast<double>(std::real(N)) - static_cast<double>(nele) < d0){
/*
  Too few electrons. Increase the chemical potential
*/

          b_ll = std::real(lambda) ;
          lambda = static_cast<z>((b_ul + b_ll)/d2) ;
        } else {
          b_ul = std::real(lambda) ;
          lambda = static_cast<z>((b_ul + b_ll)/d2) ;
          }
        }
      H = h + W + lvlshift*R + lambda*mu ;
      }

    k = c.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis).conjugate()*c.block( 0, 2*nbasis, 2*nbasis, 2*nbasis).transpose() ;

/*
  Now compare densities to check convergence.
*/
    t = p_prev - p ;
    energy = t.norm() ;
    std::cout << " rms difference in rho: " << energy << std::endl ;
    t = k_prev - k ;
    energy += t.norm() ;
    std::cout << " rms difference in kappa: " << t.norm() << std::endl ;
/*
  Check that the density has converged
*/
    if ( std::real(energy) <= thresh ) {
        std::cout << " Converged on the density " << std::endl ;
        break ;
      }

    p_prev = p ;
    k_prev = k ;
    R.block( 0, 0, 2*nbasis, 2*nbasis) = p ;
    R.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis) = I - p.conjugate() ;
    R.block( 0, 2*nbasis, 2*nbasis, 2*nbasis) = k ;
    R.block( 2*nbasis, 0, 2*nbasis, 2*nbasis) = -k.conjugate() ;
    X->contract( p, k) ;
    W = X->getG() ;
  }

/*
  Save the eigenvalues and vectors
*/

  t = (static_cast<z>(d2)*h.block( 0, 0, 2*nbasis, 2*nbasis) + W.block( 0, 0, 2*nbasis, 2*nbasis))*p ;
  energy = t.trace() ;
  std::cout << " HF Energy " << std::real(energy)/d2 << std::endl ;
  t = k.transpose()*W.block( 0, 2*nbasis, 2*nbasis, 2*nbasis) ;
  std::cout << " Pairing Energy " << std::real(t.trace())/d2 << std::endl ;
  energy += t.trace() ;
  std::cout << " Total Electronic Energy " << std::real(energy)/d2 << std::endl ;
  eig = H_diag.eigenvalues() ;
  c = H_diag.eigenvectors() ;

/*
  Clean up the memory
*/

  t.resize( 0, 0) ;
  mu.resize( 0, 0) ;
  W.resize( 0, 0) ;
  D.resize( 0, 0) ;
  G.resize( 0, 0) ;
  I.resize( 0, 0) ;
  R.resize( 0, 0) ;
  k_prev.resize( 0, 0) ;
  p_prev.resize( 0, 0) ;
  H.resize( 0, 0) ;

  ghfbdia_time.end() ;

  return std::real(energy)/d2 ;

} ;

template double ghfbdia( const Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*, const int&, const int&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::Ref<Eigen::VectorXd>, double&, double&, const int&, const int&, const double&) ;

template double ghfbdia( const Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*, const int&, const int&, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::MatrixXcd&, Eigen::Ref<Eigen::VectorXd>, cd&, cd&, const int&, const int&, const double&) ;

