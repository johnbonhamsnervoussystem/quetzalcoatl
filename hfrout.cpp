#include <cmath>
#include <complex>
#include "common.h"
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "evalm.h"
#include "hfrout.h"
#include <iostream>
#include "qtzio.h"
#include "wfn.h"
#include "solver.h"
#include <string>
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include <vector>

/* This set of routines finds a wavefunction by repeated diagonalization of the fock matrix for 
 * all the flavors of HF.  
 *
 * */
void scf_drv( common& com, std::vector<tei>& intarr, int opt) {
  /* 
     Adding HFB :
      1 - real restricted
      2 - complex restricted
      3 - real unrestricted
      4 - complex unrestricted
      5 - real generalized
      6 - complex generalized
     10 - Slater Determinant
     20 - Hartree-Fock-Bogoliubov
  */
  if ( opt / 10 == 2) {
    if ( (opt % 10) % 2 == 1 ) {
      /* Do real scf */
      real_HFB( com, intarr, opt % 10) ;
    } else if ( (opt % 10) % 2 == 0 ) {
      /* Do Complex scf */
      std::cout << " Complex HFB not yet implemented " << std::endl ;
      std::cout << " Exiting Quetzalcoatl " << std::endl ;
      exit(EXIT_FAILURE) ;
      }
  } else if ( opt / 10 == 1) {
    if ( (opt % 10) % 2 == 1 ) {
    /* Do real scf */
      real_SlaDet( com, intarr, opt % 10) ;
    } else if ( (opt % 10) % 2 == 0 ) {
    /* Do complex scf */
      cplx_SlaDet( com, intarr, opt % 10) ;
      }
    } else {
      std::cout << " Unrecognized option in scf_drv " << std::endl ;
      std::cout << " Exiting Quetzalcoatl " << std::endl ;
      exit(EXIT_FAILURE) ;
      }
  } ;

void real_HFB( common& com, std::vector<tei>& intarr, int opt) {
    /* Do the general HFB first. */
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int maxit = com.mxscfit() ;
    double thresh = com.scfthresh() ;
    double e_scf ;
    Eigen::MatrixXd h ;
    Eigen::MatrixXd s ;
    Eigen::MatrixXd moc ;
    Eigen::VectorXd eig ;
    time_dbg real_HFB_time = time_dbg("real_HFB") ;
    if ( (opt+1)/2 == 1 ) {
      std::cout << "Not yet implemented" << std::endl ;
      std::cout << "Exiting Quetzalcoatl" << std::endl ;
      exit(EXIT_FAILURE) ;
    } else if ( (opt+1)/2 == 2 ) {
      std::cout << "Not yet implemented" << std::endl ;
      std::cout << "Exiting Quetzalcoatl" << std::endl ;
      exit(EXIT_FAILURE) ;
    } else if ( (opt+1)/2 == 3 ) {
      h.resize( 4*nbas, 4*nbas) ;
      s.resize( 4*nbas, 4*nbas) ;
      moc.resize( 4*nbas, 4*nbas) ;
      eig.resize( 4*nbas) ;
      h.setZero() ;
      h.block( 0, 0, nbas, nbas) = com.getH() ;
      h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
      h.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = -h.block( 0, 0, 2*nbas, 2*nbas) ;
      s.setZero() ;
      s.block( 0, 0, nbas, nbas) = com.getS() ;
      s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
      s.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = s.block( 0, 0, 2*nbas, 2*nbas) ;
      e_scf = rghfbdia( h, s, intarr, nbas, nele, moc, eig, maxit, thresh) ;
    } else {
      std::cout << "Unrecognized option in real_HFB" << std::endl ;
      std::cout << "Exiting Quetzalcoatl" << std::endl ;
      exit(EXIT_FAILURE) ;
      }

    real_HFB_time.end() ;

    return ;

} ;

void real_SlaDet( common& com, std::vector<tei>& intarr, int opt){
/*  
    Do scf for real wavefunctions.  
      The dimension is done here but inside the routine it takes 
      care of other speficic details of the model.
    rhfdia - restricted hartree fock -> 1
    uhfdia - unrestricted hartree fock -> 2
    ghfdia - generalized hartree fock -> 3 */
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int nalp = com.nalp() ;
    int nbet = com.nbet() ;
    int maxit = com.mxscfit() ;
    double thresh = com.scfthresh() ;
    Eigen::MatrixXd h ;
    Eigen::MatrixXd s ;
    Eigen::MatrixXd ca ;
    Eigen::MatrixXd cb ;
    Eigen::VectorXd eig ;
    wfn< double, Eigen::Dynamic, Eigen::Dynamic> w ;
    time_dbg real_scf_time = time_dbg("real_scf") ;

    if ( (opt+1)/2 == 1 ) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      h = com.getH() ;
      s = com.getS() ;
      w.moc.resize( nbas, nbas) ;
      w.moc.setZero() ;
      w.eig.resize( nbas) ;
      w.e_scf = rrhfdia( h, s, intarr, nbas, nele, w.moc, w.eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << w.moc << std::endl ;
    } else if ( (opt+1)/2 == 2 ) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      w.moc.resize( nbas, 2*nbas) ;
      w.moc.setZero() ;
      h = com.getH() ;
      s = com.getS() ;
      w.eig.resize( 2*nbas) ;
      w.e_scf = ruhfdia( h, s, intarr, nbas, nalp, nbet, w.moc.block( 0, 0, nbas, nbas), w.moc.block( 0, nbas, nbas, nbas), w.eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "Alpha MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig.head(nbas) << std::endl ;
      std::cout << "Beta MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig.tail(nbas) << std::endl ;
      std::cout << "Alpha MO coefficients : " << std::endl << std::endl ;
      std::cout <<  w.moc.block( 0, 0, nbas, nbas) << std::endl ;
      std::cout << "Beta MO coefficients : " << std::endl << std::endl ;
      std::cout <<  w.moc.block( 0, nbas, nbas, nbas) << std::endl ;
    } else if ( (opt+1)/2 == 3 ) {
      h.resize( 2*nbas, 2*nbas) ;
      h.setZero() ;
      s.resize( 2*nbas, 2*nbas) ;
      s.setZero() ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      w.moc.setZero() ;
      h.block( 0, 0, nbas, nbas) = com.getH() ;
      h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
      s.block( 0, 0, nbas, nbas) = com.getS() ;
      s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
      w.eig.resize( 2*nbas) ;
      w.e_scf = rghfdia( h, s, intarr, nbas, nele, w.moc, w.eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << w.moc << std::endl ;
    } else {
      std::cout << "Unrecognized option in real_SlaDet" << std::endl ;
      std::cout << "Exiting Quetzalcoatl" << std::endl ;
      exit(EXIT_FAILURE) ;
      }

    save_slater_det(w) ;
    w.eig.resize( 0) ;
    w.moc.resize( 0, 0) ;
    s.resize( 0, 0) ;
    h.resize( 0, 0) ;

    real_scf_time.end() ;

    return ;

  } ;

void cplx_SlaDet( common& com, std::vector<tei>& intarr, int opt){
  /*Driver routine for solving self consistent wavefunctions and various things.
   First I will handle Slater-Determinants.e
    crhfdia - complex restricted hartree fock -> 1
    cuhfdia - complex unrestricted hartree fock -> 2
    cghfdia - complex generalized hartree fock -> 3 */
    int nbas = com.nbas() ;
    int nele = com.nele() ;
    int nalp = com.nalp() ;
    int nbet = com.nbet() ;
    int maxit = com.mxscfit() ;
    double thresh = com.scfthresh() ;
    double energy ;
    Eigen::MatrixXcd h ;
    Eigen::MatrixXcd s ;
    Eigen::MatrixXcd ca ;
    Eigen::MatrixXcd cb ;
    Eigen::VectorXd eig ;
    wfn< cd, Eigen::Dynamic, Eigen::Dynamic> w ;
    time_dbg cplx_scf_time = time_dbg("cplx_scf") ;

    if ( opt/2 == 1 ) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      h.setZero() ;
      s.setZero() ;
      h.real() = com.getH() ;
      s.real() = com.getS() ;
      w.moc.resize( nbas, nbas) ;
      w.moc.setZero() ;
      w.eig.resize( nbas) ;
      w.e_scf = crhfdia( h, s, intarr, nbas, nele, w.moc, w.eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << w.moc << std::endl ;
    } else if ( opt/2 == 2 ) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      h.setZero() ;
      s.setZero() ;
      h.real() = com.getH() ;
      s.real() = com.getS() ;
      w.moc.resize( nbas, 2*nbas) ;
      w.moc.setZero() ;
      w.eig.resize( 2*nbas) ;
      w.e_scf = cuhfdia( h, s, intarr, nbas, nalp, nbet, w.moc.block( 0, 0, nbas, nbas), w.moc.block( 0, nbas, nbas, nbas), w.eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "Alpha MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig.head(nbas) << std::endl ;
      std::cout << "Beta MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig.tail(nbas) << std::endl ;
      std::cout << "Alpha MO coefficients : " << std::endl << std::endl ;
      std::cout <<  w.moc.block( 0, 0, nbas, nbas) << std::endl ;
      std::cout << "Beta MO coefficients : " << std::endl << std::endl ;
      std::cout <<  w.moc.block( 0, nbas, nbas, nbas) << std::endl ;
    } else if ( opt/2 == 3 ) {
      h.resize( 2*nbas, 2*nbas) ;
      h.setZero() ;
      s.resize( 2*nbas, 2*nbas) ;
      s.setZero() ;
      w.moc.resize( 2*nbas, 2*nbas) ;
      w.moc.setZero() ;
      h.block( 0, 0, nbas, nbas) = com.getH() ;
      h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
      s.block( 0, 0, nbas, nbas) = com.getS() ;
      s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
      w.eig.resize( 2*nbas) ;
      w.e_scf = cghfdia( h, s, intarr, nbas, nele, w.moc, w.eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << w.moc << std::endl ;
    } else {
      std::cout << "Unrecognized option in cplx_SlaDet" << std::endl ;
      std::cout << "Exiting Quetzalcoatl" << std::endl ;
      exit(EXIT_FAILURE) ;
      }

    save_slater_det(w) ;
    w.eig.resize( 0) ;
    w.moc.resize( 0, 0) ;
    s.resize( 0, 0) ;
    h.resize( 0, 0) ;
    cplx_scf_time.end() ;

    return ;

  } ;

double rrhfdia( Eigen::Ref<Eigen::MatrixXd> h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, int nbasis, int nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh){

  /* Real restricted Hartree-Fock solved by repeated diagonalization. */
  Eigen::MatrixXd f ;
  Eigen::MatrixXd g ;
  Eigen::MatrixXd p ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> f_diag ;
  int iter=0 ;
  int occ ;
  double energy ;
  double ene_p=d0 ;
  double e_dif=1e0 ;
  time_dbg rrhfdia_time = time_dbg("rrhfdia") ;

  occ = nele/2 ;
  f.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p.resize( nbasis, nbasis) ;

  /* If c has something in it use it as the initial guess. */
  if( c.isZero(0) ) {
    f = h ;
  } else {
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
  } 

  while ( iter < maxit ) {
    iter += 1 ;
    f_diag.compute( f, s) ;
    c = f_diag.eigenvectors().real() ;
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
    g = p*( h + f) ;

    energy = g.trace() ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues() ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  rrhfdia_time.end() ;

  return energy ;

} ;

double crhfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh){

  /* Compelx restricted Hartree-Fock solved by repeated diagonalization. */
  typedef std::complex<double> cf ;
  Eigen::MatrixXcd f ;
  Eigen::MatrixXcd g ;
  Eigen::MatrixXcd p ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> f_diag ;
  int iter=0 ;
  int occ ;
  double energy=d0 ;
  std::complex<double> t_f ;
  double ene_p=d0 ;
  double e_dif=1e0 ;
  time_dbg crhfdia_time = time_dbg("crhfdia") ;

  occ = nele/2 ;
  f.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p.resize( nbasis, nbasis) ;

  /* If something is saved in c do an inital guess. */
  if( c.isZero(0) ) {
    f = h ;
  } else {
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
  } 

  while ( iter < maxit ) {
    iter += 1 ;
    f_diag.compute( f, s) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    ctr2er( intarr, p, g, nbasis) ;
    f = h + g ;
    g = p*(h + f) ;
    t_f = g.trace() ;
    energy = t_f.real() ;
    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
  }

  eig = f_diag.eigenvalues().real() ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  crhfdia_time.end() ;

  return energy ;

} ;

double ruhfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXd> c_a, Eigen::Ref<Eigen::MatrixXd> c_b, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh){

  /* Real unrestricted Hartree-Fock solved by repeated diagonalization. */
  Eigen::MatrixXd f_a ;
  Eigen::MatrixXd f_b ;
  Eigen::MatrixXd g ;
  Eigen::MatrixXd p_a ;
  Eigen::MatrixXd p_b ;
  Eigen::MatrixXd p_t ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> f_diag ;

  int iter=0 ;
  int nbas ;
  double energy=d0 ;
  double t_f=d0 ;
  double ene_p=d0 ;
  double e_dif=1e0 ;
  time_dbg ruhfdia_time = time_dbg("ruhfdia") ;

  nbas = 2*nbasis ;
  f_a.resize( nbasis, nbasis) ;
  f_b.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p_a.resize( nbasis, nbasis) ;
  p_b.resize( nbasis, nbasis) ;
  p_t.resize( nbasis, nbasis) ;

  /* If c_a and c_b have values use them as the initial guess */
  if ( c_a.isZero(0) ) {
    f_a = h ;
    f_b = h ;
  } else {
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
    } 

  while ( iter < maxit ) {
    iter += 1 ;
    f_diag.compute( f_a, s) ;
    c_a = f_diag.eigenvectors().real() ;
    f_diag.compute( f_b, s) ;
    c_b = f_diag.eigenvectors().real() ;
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    g = p_a*f_a ;
    t_f = g.trace() ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
    g = p_b*f_b ;
    t_f += g.trace() ;
    g = p_t*h ;
    t_f += g.trace() ;
    energy = t_f/2.0 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
    }

  f_diag.compute( f_a, s) ;
  eig.head(nbasis) = f_diag.eigenvalues() ;
  f_diag.compute( f_b, s) ;
  eig.tail(nbasis) = f_diag.eigenvalues() ;
  p_t.resize( 0, 0) ;
  p_b.resize( 0, 0) ;
  p_a.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f_b.resize( 0, 0) ;
  f_a.resize( 0, 0) ;

  ruhfdia_time.end() ;

  return energy ;

} ;

double cuhfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXcd> c_a, Eigen::Ref<Eigen::MatrixXcd> c_b, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh){

  /* Complex unrestricted Hartree-Fock solved by repeated diagonalization. */
  typedef std::complex<double> cf ;
  Eigen::MatrixXcd f_a ;
  Eigen::MatrixXcd f_b ;
  Eigen::MatrixXcd g ;
  Eigen::MatrixXcd p_a ;
  Eigen::MatrixXcd p_b ;
  Eigen::MatrixXcd p_t ;
  Eigen::MatrixXcd temp ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> f_diag ;
  int iter=0 ;
  int nbas ;
  double energy=d0 ;
  std::complex<double> t_f ;
  double ene_p=d0 ;
  double e_dif=1e0 ;
  time_dbg cuhfdia_time = time_dbg("cuhfdia") ;

  nbas = 2*nbasis ;
  f_a.resize( nbasis, nbasis) ;
  f_b.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p_a.resize( nbasis, nbasis) ;
  p_b.resize( nbasis, nbasis) ;
  p_t.resize( nbasis, nbasis) ;
  temp.resize( nbasis, nbasis) ;

  /* If c_a and c_b are filled use them for an inital guess */
  if ( c_a.isZero(0) ) {
    f_a = h ;
    f_b = h ;
  } else {
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
    } 

  while ( iter < maxit ) {
    iter += 1 ;
    f_diag.compute( f_a, s) ;
    c_a = f_diag.eigenvectors() ;
    f_diag.compute( f_b, s) ;
    c_b = f_diag.eigenvectors() ;
    p_a = c_a.block( 0, 0, nbasis, nalp)*c_a.block( 0, 0, nbasis, nalp).adjoint() ;
    p_b = c_b.block( 0, 0, nbasis, nbet)*c_b.block( 0, 0, nbasis, nbet).adjoint() ;
    p_t = p_a + p_b ;
    ctr2eu( intarr, p_t, p_a, g, nbasis) ;
    f_a = h + g ;
    g = p_a*f_a ;
    temp = g ;
    ctr2eu( intarr, p_t, p_b, g, nbasis) ;
    f_b = h + g ;
    g = p_b*f_b ;
    temp = temp + g ;
    g = p_t*h ;
    temp = temp + g ;
    t_f =  temp.trace() ;
    energy = t_f.real()/2.0 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 5 && e_dif < thresh ) { break ;}
    }

  f_diag.compute( f_a, s) ;
  eig.head(nbasis) = f_diag.eigenvalues().real() ;
  f_diag.compute( f_b, s) ;
  eig.tail(nbasis) = f_diag.eigenvalues().real() ;
  temp.resize( 0, 0) ;
  p_t.resize( 0, 0) ;
  p_b.resize( 0, 0) ;
  p_a.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f_b.resize( 0, 0) ;
  f_a.resize( 0, 0) ;

  cuhfdia_time.end() ;

  return energy ;

} ;

double rghfbdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh) {

/* 
  There are two convergence criteria here to achieve self-consistency. 

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
    maxit - maxiterations before convergence is stopped
    thresh - threshold for convergence o the density

  Local :
    H - HFB Hamiltonian
    p - Stores the normal density
    k - Stores the abnormal density
    W - This is a container for the lowest eigenvectors of the HFB quasiparticles
    R - generalized density matrix
    G - Self-Consistent Field 
    D - Pairing Field 
    C - Container for Self-Consisten and Pairing fields
    mu - anti-symmetric overlap for updating chemical potential
    iter_d - iterations on density convergence
    iter_N - iterations on particle number
    lambda - chemical potential
    update_lambda - magnitude by which to change the chemical potential
    H_diag - Eigensolver object
    Ro - orthogonal density
    t - scratch space
    energy - final HFB energy/ this also stores the rms error of the density
	in order to check convergence
    N - number of particles
    N_p - number of particles from a previous iteration

*/

  Eigen::MatrixXd H ;
  Eigen::MatrixXd p ;
  Eigen::MatrixXd k ;
  Eigen::MatrixXd W ;
  Eigen::MatrixXd R ;
  Eigen::MatrixXd G ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXd C ;
  Eigen::MatrixXd mu ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> H_diag ;
  Eigen::MatrixXd Ro ;
  Eigen::MatrixXd t ;
/*
*/
  int iter_N = 0 ;
  int iter_d = 0 ;
  double lambda=-d1 ;
  double update_lambda=d1 ;
  double energy=d0 ;
  double N ;
  double N_p=-N ;
  double conv=d0 ;
  time_dbg rghfbdia_time = time_dbg("rghfbdia") ;

  H.resize( nbasis*4, nbasis*4) ;
  p.resize( nbasis*2, nbasis*2) ;
  k.resize( nbasis*2, nbasis*2) ;
  W.resize( nbasis*4, nbasis*2) ;
  R.resize( nbasis*4, nbasis*4) ;
  G.resize( nbasis*2, nbasis*2) ;
  D.resize( nbasis*2, nbasis*2) ;
  mu.resize( nbasis*4, nbasis*4) ;
  Ro.resize( nbasis*4, nbasis*4) ;

  C.resize( nbasis*4, nbasis*4) ;
  t.resize( nbasis*2, nbasis*2) ;

  /* Initial guess for rho and kappa */
  if ( c.isZero(0) ) {
    H = h + C.setRandom() ;
    p.setZero() ;
    k.setZero() ;
  } else {
    W = c.block(0, 0, 4*nbasis, 2*nbasis) ;
    R = W*W.adjoint() ;
    p = R.block(0, 0, 2*nbasis, 2*nbasis) ;
    k = R.block(0, 2*nbasis, 2*nbasis, 2*nbasis) ;
    ctr2eg( intarr, p, G, nbasis) ;
    ctrPairg( intarr, k, D, nbasis) ;
    C.setZero() ;
    C.block( 0, 0, 2*nbasis, 2*nbasis) = G ;
    C.block( 0, 2*nbasis, 2*nbasis, 2*nbasis) = D/d2 ;
    C.block( 2*nbasis, 0, 2*nbasis, 2*nbasis) = -(D.conjugate())/d2 ;
    C.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis) = -G.conjugate() ;
    H.setZero() ;
    H = h + C + lambda*mu ;
    }

  mu.setZero() ;
  mu.block( 0, 0, 2*nbasis, 2*nbasis) = -s.block( 0, 0, 2*nbasis, 2*nbasis) ;
  mu.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis) = s.block( 0, 0, 2*nbasis, 2*nbasis) ;

  while ( iter_d < maxit ) {
    iter_d += 1 ;
    std::cout << std::endl << "Self-Consistency iteration number: " << iter_d << std::endl ;
    iter_N = 0 ;
    while ( iter_N < maxit ){
      iter_N += 1 ;
      std::cout << "  Particle Number Iteration: " << iter_N << std::endl ;
      std::cout << "    chemical potential: " << lambda << std::endl ;
      H_diag.compute( H, s) ;
      c = H_diag.eigenvectors().real() ;
      W = c.block(0, 0, 4*nbasis, 2*nbasis) ;
      R = W*W.adjoint() ;
      Ro = R*s ;
      N = Ro.block( 0, 0, 2*nbasis, 2*nbasis).trace() ;
      std::cout << "    Particle Number: " << N << std::endl ;
      if  ( abs(N - static_cast<double>(nele)) < 0.1 ){
        break ;
      } else if ( N > static_cast<double>(nele)){
        lambda += -update_lambda ;
      } else {
        lambda += update_lambda ;
      }
      H = h + C + lambda*mu ;
    }
/*
  Now compare densities to check convergence.
*/
    t = p - R.block(0, 0, 2*nbasis, 2*nbasis) ;
    energy = sqrt(t.squaredNorm()) ;
    t = k - R.block(0, 2*nbasis, 2*nbasis, 2*nbasis) ;
    energy += sqrt(t.squaredNorm()) ;
    std::cout << "  rms difference in the densities: " << energy << std::endl ;

/*  Update the normal and abnormal density before we check convergence. */

    p = R.block(0, 0, 2*nbasis, 2*nbasis) ;
    k = R.block(0, 2*nbasis, 2*nbasis, 2*nbasis) ;

/*
  It's not immediately clear to me how to optimize the chemical potential.
  For now I will have it make steps of size 1.  
*/

/* Check that the density has converged */
    if ( energy <= thresh){
/*    Check that the particle number has converged */
      if ( abs(N - N_p) < 0.1 ) {
	break ;
      } else {
        N_p = N ;
        }
      }
/*
  If we have not converged then build a new Hamiltonian
*/
    ctr2eg( intarr, p, G, nbasis) ;
    ctrPairg( intarr, k, D, nbasis) ;
    C.setZero() ;
    C.block( 0, 0, 2*nbasis, 2*nbasis) = G ;
    C.block( 0, 2*nbasis, 2*nbasis, 2*nbasis) = D/d2 ;
    C.block( 2*nbasis, 0, 2*nbasis, 2*nbasis) = -(D.conjugate())/d2 ;
    C.block( 2*nbasis, 2*nbasis, 2*nbasis, 2*nbasis) = -G.conjugate() ;
    H.setZero() ;
    H = h + C + lambda*mu ;
  }

  t = (h.block( 0, 0, 2*nbasis, 2*nbasis) + G/d2 - lambda*s.block( 0, 0, 2*nbasis, 2*nbasis))*p ;
  energy = t.trace() ;
  t = k.conjugate()*D ;
  energy += t.trace()/4.0e0 ;
  std::cout << " HFB energy: " << energy << std::endl ;

  t.resize( 0, 0) ;
  Ro.resize( 0, 0) ;
  mu.resize( 0, 0) ;
  D.resize( 0, 0) ;
  G.resize( 0, 0) ;
  R.resize( 0, 0) ;
  W.resize( 0, 0) ;
  k.resize( 0, 0) ;
  p.resize( 0, 0) ;
  H.resize( 0, 0) ;

  rghfbdia_time.end() ;

  return energy ;

} ;

double rghfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh) {

  /* Real generalized Hartree-Fock solved by repeated diagonalization. */
  Eigen::MatrixXd f ;
  Eigen::MatrixXd g ;
  Eigen::MatrixXd p ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> f_diag ;
  int iter=0 ;
  int nbas ;
  double energy=d0 ;
  double ene_p=d0 ;
  double e_dif=1e0 ;
  time_dbg rghfdia_time = time_dbg("rghfdia") ;

  nbas = nbasis*2 ;
  f.resize( nbas, nbas) ;
  g.resize( nbas, nbas) ;
  p.resize( nbas, nbas) ;

  /* If something is stored in c use it for an inital guess */
  if ( c.isZero(0) ) {
    f = h ;
  } else {
    c = f_diag.eigenvectors().real() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h + g ;
  } 

  while ( iter < maxit ) {
    iter += 1 ;
    f_diag.compute( f, s) ;
    c = f_diag.eigenvectors().real() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h + g ;
    g = p*( h + f) ;
    energy = g.trace()/d2 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 3 && e_dif < thresh ) { break ;}
    }

  std::cout << " energy = " << energy << std::endl ;

  eig = f_diag.eigenvalues() ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  rghfdia_time.end() ;

  return energy ;

} ;

double cghfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh){

  /* Compelx restricted Hartree-Fock solved by repeated diagonalization.
     This may be eventually reduces to a single block diagonalization
     but for now the electron contraction routines are handled by differnt
     algorithms.  */
  Eigen::MatrixXcd f ;
  Eigen::MatrixXcd g ;
  Eigen::MatrixXcd p ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> f_diag ;
  int iter=0 ;
  int nbas ;
  double energy=d0 ;
  cd t_f ;
  double ene_p=d0 ;
  double e_dif=1e0 ;
  time_dbg cghfdia_time = time_dbg("cghfdia") ;

  nbas = nbasis*2 ;
  f.resize( nbas, nbas) ;
  g.resize( nbas, nbas) ;
  p.resize( nbas, nbas) ;

  if ( c.isZero(0) ) {
    f = h ;
  } else {
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h + g ;
  } 

  while ( iter < maxit ) {
    iter += 1 ;
    f_diag.compute( f, s) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbas, nele)*c.block( 0, 0, nbas, nele).adjoint() ;
    ctr2eg( intarr, p, g, nbasis) ;
    f = h + g ;
    g = p*( h + f) ;
    t_f = g.trace() ;
    energy = t_f.real()/d2 ;

    e_dif = std::abs(ene_p - energy) ;
    ene_p = energy ;
    if ( iter > 3 && e_dif < thresh ) { break ;}
    }

  eig = f_diag.eigenvalues().real() ;
  p.resize( 0, 0) ;
  g.resize( 0, 0) ;
  f.resize( 0, 0) ;

  cghfdia_time.end() ;

  return energy ;

} ;

