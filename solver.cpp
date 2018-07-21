#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "binio.h"
#include "common.h"
#include "hfwfn.h"
#include "solver.h"
#include "tei.h"
#include "util.h"
#include "evalm.h"

/*
 * Notes for Solver:
 *
 *   - KT 06/2018 Added ahm_exp here because it makes use of the 
 *     Eigenvalues package. 
 *   - Find a real convergence criteria for ahm_exp based on repreated
 *     mutliplication
 *
 * */

void canort( Eigen::Ref<Eigen::MatrixXd> const s, Eigen::Ref<Eigen::MatrixXcd> xs , int& dim) {
  /*  
 *    Find an orthogonal transformation using canonical orthogonalization for a real matrix
 *    */

  int indx=-1 ;
  double thresh = 1e-7 ;
  std::complex<double> eig ;
  Eigen::EigenSolver<Eigen::MatrixXd> s_diag( s, true) ;
  Eigen::MatrixXcd s_U ;
  Eigen::VectorXcd s_eig ;

  /* Allocate some space to build our transformation matrix. */

  s_eig.resize( dim) ;
  s_U.resize( dim, dim) ;

  s_eig = s_diag.eigenvalues() ;

  s_U = s_diag.eigenvectors() ;

  xs = s_U.householderQr().householderQ() ;

  for (int i = 0; i < dim; i++ ) {
    eig = s_eig(i) ;
    if ( std::real(eig) <= thresh && std::imag(eig) <= thresh ) {
      std::cout << "Linear dependency found in canort.  Reducing basis." << std::endl ;
      std::cout << s_eig(i) << " is less than the threshhold" << std::endl ;
      std::cout << "Reducing basis" << std::endl ;
    } else {
      indx ++ ;
    }
    xs.col(indx) = xs.col(i)/std::sqrt(eig) ;
  }
 
  if ( indx != dim-1 ){
    xs.resize(dim,indx) ;
  }

  dim = indx + 1 ;
 
  /* Deallocate  */

  s_eig.resize( 0) ;
  s_U.resize( 0, 0) ;

  return ;

} ;

void canort( Eigen::Ref<Eigen::MatrixXcd> const s, Eigen::Ref<Eigen::MatrixXcd> xs , int& dim) {
  /*  
 *    Find an orthogonal transformation using canonical orthogonalization for a complex matrix
 *    */

  int indx=-1 ;
  double thresh = 1e-7 ;
  std::complex<double> eig ;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> s_diag( s, true) ;
  Eigen::MatrixXcd s_U ;
  Eigen::VectorXcd s_eig ;

  /* Allocate some space to build our transformation matrix. */

  s_eig.resize( dim) ;
  s_U.resize( dim, dim) ;

  s_eig = s_diag.eigenvalues() ;

  s_U = s_diag.eigenvectors() ;

  xs = s_U.householderQr().householderQ() ;

  for (int i = 0; i < dim; i++ ) {
    eig = s_eig(i) ;
    if ( abs(std::real(eig)) <= thresh && abs(std::imag(eig)) <= thresh ) {
      std::cout << "Linear dependency found in canort.  Reducing basis." << std::endl ;
      std::cout << s_eig(i) << " is less than the threshhold" << std::endl ;
      std::cout << "Reducing basis" << std::endl ;
    } else {
      indx ++ ;
    }
    xs.col(indx) = xs.col(i)/std::sqrt(eig) ;
  }
 
//  if ( indx != dim-1 ){
//    xs.resize(4,indx) ;
//  }

  dim = indx + 1 ;
 
  /* Deallocate  */

  s_eig.resize( 0) ;
  s_U.resize( 0, 0) ;

  return ;

} ;

void ahm_exp( Eigen::Ref<Eigen::MatrixXcd> x, Eigen::Ref<Eigen::MatrixXcd> u, int dim,  int opt) {
/*
 *  Anti-Hermetian Matrix EXPonential
 *
 *  Given an anti-hermetian matrix in an exponential, build a unitary
 *  matrix. 
 *
 *  This is done two ways. The first is by diagonalizing the 
 *  anti-hermetian matrix
 *
 *     X = iVdVt     VtV = 1
 *
 *  U = exp(X) = exp(iVdVt) = V*exp(id)*Vt
 *
 *  The other method is by the definition of a matrix exponential
 *
 *  U = exp(X) = 1 +  sum_{n=1}^{infty} X^{n}/n!
 *
 *  where 1 is the unit matrix.
 *
 *
 * */ 
  int n = 0 ;
  double denom=1.0 ;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> xs ;
  Eigen::MatrixXcd d ;
  Eigen::MatrixXcd tmp ;

  tmp.resize( dim, dim) ;
  d.resize( dim, dim) ;
  d.setIdentity() ;

  if ( opt == 0 ) {
    /* Default behavior.  Get the eigenvalues and vectors of x */
    xs.compute( x) ;

    /* Build d */
    for ( int i = 0; i < dim; i++ ){
      d( i, i) = std::exp(xs.eigenvalues()[i]) ;
    }

    /* Build U by Vtexp(idelta)V */
    tmp = d*xs.eigenvectors() ;
    u = xs.eigenvectors().adjoint()*tmp ;

  } else if ( opt == 1 ) {
    /* Find U by repeated multiplication */
      u.setZero() ;

      /* What should the convergence criteria be? */
      while( n < 15 ) {
        u += d ;
        n ++ ;
        tmp = d*x ;
        denom = 1.0/static_cast<double>(n) ;
        d = denom*tmp ;
      }

  } else {

    std::cout << "Unrecgonized opt in ahm_exp" << std::endl ;
    std::cout << "Returning Unit matrix" << std::endl ;
    u.setIdentity() ;

  }

  return ;

} ;

void trci( common& com, std::vector<hfwfn>& det, Eigen::Ref<Eigen::MatrixXcd> H, std::vector<tei>& intarr) {

/*
 * This routine accepts N determinants into the vector. These determinants
 * are expanded in the "Time Reversal" Basis of Complex Conjugation(K),
 * Spin Flip(F), and Time Reversal(T).  This will generate 4N
 * determinants.
 * */

  Eigen::MatrixXcd CI_s ;
  Eigen::MatrixXcd CI_h ;
  Eigen::MatrixXcd tmo ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ci ;
  hfwfn bra ;
  hfwfn ket ;
  int CI_d ;
  int idt_lb ;
  int jdt_lb ;
  int idt_ty ;
  int jdt_ty ;

  CI_d = 4*det.size() ;

  CI_s.resize( CI_d, CI_d) ;
  CI_h.resize( CI_d, CI_d) ;

  /* Build a hfwfn container for calculations */
  tmo.resize( 2*com.nbas(), 2*com.nbas()) ;
  tmo.setIdentity() ;
  bra.fil_mos( com.nbas(), tmo, 6) ;
  ket.fil_mos( com.nbas(), tmo, 6) ;

  for ( int i=0; i < CI_d; i++ ) {
    idt_lb = i/4 ;
    idt_ty = i % 4 ;

    if ( idt_ty == 0 ){
    /* |Phi> */
      det[idt_lb].get_mos( tmo) ;
      bra.set_mos( tmo) ;

    } else if ( idt_ty == 1 ){
    /* K|Phi> */
      K_op( det[idt_lb], tmo, com.nbas()) ;
      bra.set_mos( tmo) ;

    } else if ( idt_ty == 2 ){
    /* F|Phi> */
      F_op( det[idt_lb], tmo, com.nbas()) ;
      bra.set_mos( tmo) ;

    } else if ( idt_ty == 3 ){
    /* T|Phi> */
      T_op( det[idt_lb], tmo, com.nbas()) ;
      bra.set_mos( tmo) ;

    }

    for ( int j=0; j < CI_d; j++ ) {
      jdt_lb = j/4 ;
      jdt_ty = j % 4 ;

      if ( jdt_ty == 0 ){
      /* |Phi> */
        det[jdt_lb].get_mos( tmo) ;
        ket.set_mos( tmo) ;
 
      } else if ( jdt_ty == 1 ){ 
      /* K|Phi> */
        K_op( det[jdt_lb], tmo, com.nbas()) ;
        ket.set_mos( tmo) ;

      } else if ( jdt_ty == 2 ){
      /* F|Phi> */
        F_op( det[jdt_lb], tmo, com.nbas()) ;
        ket.set_mos( tmo) ;

      } else if ( jdt_ty == 3 ){
      /* T|Phi> */
        T_op( det[jdt_lb], tmo, com.nbas()) ;
        ket.set_mos( tmo) ;

      }
      /* Build <phi|H|psi> and <phi|psi> */
        CI_h( i, j) = fockop( com, H, intarr, bra, ket, CI_s( i, j)) ;
      }
    }

  std::cout << " CI_H " << std::endl << CI_h << std::endl ;
  std::cout << " CI_S " << std::endl << CI_s << std::endl ;

  ci.compute( CI_h, CI_s) ;

  std::cout << " ci.eval " << std::endl << ci.eigenvalues() << std::endl ;
  std::cout << " ci.evec " << std::endl << ci.eigenvectors() << std::endl ;

  tmo.resize( 0, 0) ;
  CI_h.resize( 0, 0) ;
  CI_s.resize( 0, 0) ;

  return ;

} ;

void trci( common& com, std::vector<hfwfn>& det, Eigen::Ref<Eigen::MatrixXcd> H, std::vector<tei>& intarr, std::string& trd, std::string& fop) {

/*
 * The last two arguments given are files for writing intermediate quantities in the derivatives.
 * */

  Eigen::MatrixXcd CI_s ;
  Eigen::MatrixXcd CI_h ;
  Eigen::MatrixXcd tmo ;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ci ;
  hfwfn bra ;
  hfwfn ket ;
  int CI_d ;
  int idt_lb ;
  int jdt_lb ;
  int idt_ty ;
  int jdt_ty ;
  std::ofstream tdfile(trd, std::ios::out | std::ios::binary ) ;
  std::ofstream fofile(fop, std::ios::out | std::ios::binary ) ;

  /* Open our files for writing to. */
  CI_d = 4*det.size() ;

  CI_s.resize( CI_d, CI_d) ;
  CI_h.resize( CI_d, CI_d) ;

  /* Build a hfwfn container for calculations */
  tmo.resize( 2*com.nbas(), 2*com.nbas()) ;
  tmo.setIdentity() ;
  bra.fil_mos( com.nbas(), tmo, 6) ;
  ket.fil_mos( com.nbas(), tmo, 6) ;

  for ( int i=0; i < CI_d; i++ ) {
    idt_lb = i/4 ;
    idt_ty = i % 4 ;

    if ( idt_ty == 0 ){
    /* |Phi> */
      det[idt_lb].get_mos( tmo) ;
      bra.set_mos( tmo) ;

    } else if ( idt_ty == 1 ){
    /* K|Phi> */
      K_op( det[idt_lb], tmo, com.nbas()) ;
      bra.set_mos( tmo) ;

    } else if ( idt_ty == 2 ){
    /* F|Phi> */
      F_op( det[idt_lb], tmo, com.nbas()) ;
      bra.set_mos( tmo) ;

    } else if ( idt_ty == 3 ){
    /* T|Phi> */
      T_op( det[idt_lb], tmo, com.nbas()) ;
      bra.set_mos( tmo) ;

    }

    for ( int j=0; j < CI_d; j++ ) {
      jdt_lb = j/4 ;
      jdt_ty = j % 4 ;

      if ( jdt_ty == 0 ){
      /* |Phi> */
        det[jdt_lb].get_mos( tmo) ;
        ket.set_mos( tmo) ;
 
      } else if ( jdt_ty == 1 ){ 
      /* K|Phi> */
        K_op( det[jdt_lb], tmo, com.nbas()) ;
        ket.set_mos( tmo) ;

      } else if ( jdt_ty == 2 ){
      /* F|Phi> */
        F_op( det[jdt_lb], tmo, com.nbas()) ;
        ket.set_mos( tmo) ;

      } else if ( jdt_ty == 3 ){
      /* T|Phi> */
        T_op( det[jdt_lb], tmo, com.nbas()) ;
        ket.set_mos( tmo) ;

      }
      /* Build <phi|H|psi> and <phi|psi> */
        CI_h( i, j) = fockop( com, H, intarr, bra, ket, CI_s( i, j), tdfile, fofile) ;
      }
    }

  std::cout << " CI_H " << std::endl << CI_h << std::endl ;
  std::cout << " CI_S " << std::endl << CI_s << std::endl ;

  ci.compute( CI_h, CI_s) ;

  std::cout << " ci.eval " << std::endl << ci.eigenvalues() << std::endl ;
  std::cout << " ci.evec " << std::endl << ci.eigenvectors() << std::endl ;

  tmo.resize( 0, 0) ;
  CI_h.resize( 0, 0) ;
  CI_s.resize( 0, 0) ;

  tdfile.close() ;
  fofile.close() ;

  return ;

} ;
