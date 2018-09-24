#include <complex>
#include "common.h"
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "evalm.h"
#include "hfrout.h"
#include <iostream>
#include "qtzio.h"
#include "sladet.h"
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
  /* Right now we just do Slater determinants. Go figure. */
  if ( (opt / 10) % 2 == 0 ){
  /* Do real scf */
    real_SlaDet( com, intarr, opt % 10) ;
  } else {
  /* Do complex scf */
    cplx_SlaDet( com, intarr, opt % 10) ;
    }
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
    double energy ;
    Eigen::MatrixXd h ;
    Eigen::MatrixXd s ;
    Eigen::MatrixXd ca ;
    Eigen::MatrixXd cb ;
    Eigen::VectorXd eig ;
    time_dbg real_scf_time = time_dbg("real_scf") ;
    sladet< double, Eigen::Dynamic, Eigen::Dynamic> w ;

    if ( opt == 1) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      w.moc.resize( nbas, nbas) ;
      w.moc.setZero() ;
      h = com.getH() ;
      s = com.getS() ;
      w.eig.resize( nbas) ;
      w.e_scf = rrhfdia( h, s, intarr, nbas, nele, w.moc, w.eig, maxit, thresh) ;
      save_slater_det( w) ;
      std::cout << "Mean Field Energy : " << w.e_scf + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << w.eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << w.moc << std::endl ;
      w.eig.resize( 0) ;
      w.moc.resize( 0, 0) ;
      s.resize( 0, 0) ;
      h.resize( 0, 0) ;
    } else if ( opt == 2) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      ca.resize( nbas, nbas) ;
      cb.resize( nbas, nbas) ;
      ca.setZero() ;
      cb = ca ;
      h = com.getH() ;
      s = com.getS() ;
      eig.resize( 2*nbas) ;
      energy = ruhfdia( h, s, intarr, nbas, nalp, nbet, ca, cb, eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << energy + com.nrep() << std::endl ;
      std::cout << "Alpha MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig.head(nbas) << std::endl ;
      std::cout << "Beta MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig.tail(nbas) << std::endl ;
      std::cout << "Alpha MO coefficients : " << std::endl << std::endl ;
      std::cout << ca << std::endl ;
      std::cout << "Beta MO coefficients : " << std::endl << std::endl ;
      std::cout << cb << std::endl ;
      eig.resize( 0) ;
      cb.resize( 0, 0) ;
      ca.resize( 0, 0) ;
      s.resize( 0, 0) ;
      h.resize( 0, 0) ;
    } else if ( opt == 3) {
      h.resize( 2*nbas, 2*nbas) ;
      h.setZero() ;
      s.resize( 2*nbas, 2*nbas) ;
      s.setZero() ;
      ca.resize( 2*nbas, 2*nbas) ;
      ca.setZero() ;
      h.block( 0, 0, nbas, nbas) = com.getH() ;
      h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
      s.block( 0, 0, nbas, nbas) = com.getS() ;
      s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
      eig.resize( 2*nbas) ;
      energy = rghfdia( h, s, intarr, nbas, nele, ca, eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << energy + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << ca << std::endl ;
      eig.resize( 0) ;
      ca.resize( 0, 0) ;
      s.resize( 0, 0) ;
      h.resize( 0, 0) ;
      }

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
    time_dbg cplx_scf_time = time_dbg("cplx_scf") ;

    if ( opt == 1) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      ca.resize( nbas, nbas) ;
      h.setZero() ;
      s.setZero() ;
      ca.setZero() ;
      h.real() = com.getH() ;
      s.real() = com.getS() ;
      eig.resize( nbas) ;
      energy = crhfdia( h, s, intarr, nbas, nele, ca, eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << energy + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << ca << std::endl ;
      eig.resize( 0) ;
      ca.resize( 0, 0) ;
      s.resize( 0, 0) ;
      h.resize( 0, 0) ;
    } else if ( opt == 2) {
      h.resize( nbas, nbas) ;
      s.resize( nbas, nbas) ;
      ca.resize( nbas, nbas) ;
      cb.resize( nbas, nbas) ;
      h.setZero() ;
      s.setZero() ;
      ca.setZero() ;
      cb = ca ;
      h.real() = com.getH() ;
      s.real() = com.getS() ;
      eig.resize( 2*nbas) ;
      energy = cuhfdia( h, s, intarr, nbas, nalp, nbet, ca, cb, eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << energy + com.nrep() << std::endl ;
      std::cout << "Alpha MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig.head(nbas) << std::endl ;
      std::cout << "Beta MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig.tail(nbas) << std::endl ;
      std::cout << "Alpha MO coefficients : " << std::endl << std::endl ;
      std::cout << ca << std::endl ;
      std::cout << "Beta MO coefficients : " << std::endl << std::endl ;
      std::cout << cb << std::endl ;
      eig.resize( 0) ;
      cb.resize( 0, 0) ;
      ca.resize( 0, 0) ;
      s.resize( 0, 0) ;
      h.resize( 0, 0) ;
    } else if ( opt == 3) {
      h.resize( 2*nbas, 2*nbas) ;
      h.setZero() ;
      s.resize( 2*nbas, 2*nbas) ;
      s.setZero() ;
      ca.resize( 2*nbas, 2*nbas) ;
      ca.setZero() ;
      h.block( 0, 0, nbas, nbas).real() = com.getH() ;
      h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
      s.block( 0, 0, nbas, nbas).real() = com.getS() ;
      s.block( nbas, nbas, nbas, nbas) = s.block( 0, 0, nbas, nbas) ;
      eig.resize( 2*nbas) ;
      energy = cghfdia( h, s, intarr, nbas, nele, ca, eig, maxit, thresh) ;
      std::cout << "Mean Field Energy : " << energy + com.nrep() << std::endl ;
      std::cout << "MO Eigenvalues : " << std::endl << std::endl ;
      std::cout << eig << std::endl ;
      std::cout << "MO coefficients : " << std::endl << std::endl ;
      std::cout << ca << std::endl ;
      eig.resize( 0) ;
      ca.resize( 0, 0) ;
      s.resize( 0, 0) ;
      h.resize( 0, 0) ;
      }

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
  time_dbg rrhfdia_time = time_dbg("rrhfdia_time") ;

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
    f = g ;
    oao ( nbasis, f, s) ;
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
  time_dbg crhfdia_time = time_dbg("crhfdia_time") ;

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
  time_dbg ruhfdia_time = time_dbg("ruhfdia_time") ;

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
  time_dbg cuhfdia_time = time_dbg("cuhfdia_time") ;

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

double rghfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& maxit, double& thresh){

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
  time_dbg rghfdia_time = time_dbg("rghfdia_time") ;

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
  time_dbg cghfdia_time = time_dbg("cghfdia_time") ;

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

