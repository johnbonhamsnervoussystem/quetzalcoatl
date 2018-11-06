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
  Driver routine for projected methods.  
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

  int nbas = com.nbas() ;
  std::vector<tei> oaoint ;
  Eigen::MatrixXd Tmp ;
  Eigen::MatrixXd Xs ;
  Eigen::MatrixXcd H ;
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;

/*
  Load the basis transformation
*/

  w.moc.resize( 4*nbas, 4*nbas) ;
  w.eig.resize( 4*nbas) ;
  Tmp.resize( nbas, nbas) ;
  H.resize( 4*nbas, 4*nbas) ;
  Xs.resize( nbas, nbas) ;

  Xs = com.getXS() ;
  H.block( 0, 0, nbas, nbas).real() = com.getH() ;
  load_wfn( w) ;

  oao( Tmp, H.block( 0, 0, nbas, nbas).real(), Xs) ;
  H.block( nbas, nbas, nbas, nbas) = H.block( 0, 0, nbas, nbas) ;
  H.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = H.block( 0, 0, 2*nbas, 2*nbas) ;
/*
  Each block of V and U is made of four components and must be 
  transformed as such.
*/
  oao( Tmp, w.moc.block( 0, 0, nbas, nbas), Xs) ;
  oao( Tmp, w.moc.block( 0, nbas, nbas, nbas), Xs) ;
  oao( Tmp, w.moc.block( nbas, 0, nbas, nbas), Xs) ;
  oao( Tmp, w.moc.block( nbas, nbas, nbas, nbas), Xs) ;
/*
*/
  oao( Tmp, w.moc.block( 2*nbas, 0, nbas, nbas), Xs) ;
  oao( Tmp, w.moc.block( 2*nbas, nbas, nbas, nbas), Xs) ;
  oao( Tmp, w.moc.block( 3*nbas, 0, nbas, nbas), Xs) ;
  oao( Tmp, w.moc.block( 3*nbas, nbas, nbas, nbas), Xs) ;

  w.moc.block( 0, 2*nbas, 2*nbas, 2*nbas) = w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas).conjugate() ;
  w.moc.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) = w.moc.block( 0, 0, 2*nbas, 2*nbas).conjugate() ;

  oao( nbas, intarr, oaoint, Xs) ;

  if( true){
     if ( false){
//       proj_HFB_real() ;
     } else {
       proj_HFB_cplx( com, H, oaoint, w) ;
       }
  } else {
    ;
//    proj_HF
    }

  return ;

} ;

void proj_HFB_cplx( common& com, Eigen::Ref<Eigen::MatrixXcd>H, std::vector<tei>& oaoint, Eigen::Ref<Eigen::MatrixXcd> w){
/*
  Projected HFB
*/
  time_dbg proj_HFB_cplx_time = time_dbg("proj_HFB_cplx") ;

  Eigen::MatrixXcd m ;
  Eigen::MatrixXcd V ;
  Eigen::MatrixXcd U ;
  Eigen::MatrixXcd Uinv ;

  /* First things first, read the converged HFB wavefunciton from another job */
  wfn < cd, Eigen::Dynamic, Eigen::Dynamic> w;

  w.moc.resize( 4*nbas, 4*nbas) ;
  m.resize( 4*nbas, 4*nbas) ;
  V.resize( 2*nbas, 2*nbas) ;
  U.resize( 2*nbas, 2*nbas) ;
  Uinv.resize( 2*nbas, 2*nbas) ;

  load_wfn( w) ;

/* 
   W should have
   
     V U^{*}
     U V^{*}
*/

  V = w.moc.block( 0, 0, 2*nbas, 2*nbas) ;
  U = w.moc.block( 2*nbas, 0, 2*nbas, 2*nbas) ;
  Uinv = U.inverse() ;

  /* Let's check our implementation of the pfaffian first */

  m.block( 2*nbas, 0, 2*nbas, 2*nbas).real() =  Eigen::MatrixXd::Identity( 2*nbas, 2*nbas)  ;
  m.block( 0, 2*nbas, 2*nbas, 2*nbas).real() = - Eigen::MatrixXd::Identity( 2*nbas, 2*nbas)  ;
  U = V*Uinv ;
  m.block( 0, 0, 2*nbas, 2*nbas) =  U ;
  m.block( 2*nbas, 2*nbas, 2*nbas, 2*nbas) =  -U.conjugate() ;

  std::cout << m << std::endl ;
  
  std::cout << pfaffian( m) << std::endl ;

  proj_HFB_cplx_time.end() ;

  return ; 

} ;

void e_phf( common& com, hfwfn& refd, Eigen::Ref<Eigen::MatrixXcd> H, std::vector<tei>& intarr) {

/* 
 * Return the PHF energy given a determinant
 *
 * E = (int dOmega x(Omega) <0|H|Omega>/(int dOmega x(Omega) )
 *
 * */


/* 
 * Developing notes: I will start with a PAV that does not require
 * coefficients for spin states or complex conjugation.
 *
 * */

 // Variables

  int nbasis ;
  double Intg = 0.0 ;
  cd w_val ;
  cd energy ;
  cd ovl ;
  cd n_k = cd( 0.0, 0.0) ;
  cd h_k = cd( 0.0, 0.0) ;
  Eigen::MatrixXcd tmp ;

 // One reference determinant. One for storage space.

  hfwfn rotd ;
 
 // integration grid.

  std::vector<double> w_a ;
  std::vector<double> w_b ;
  std::vector<double> w_g ;
  std::vector<double> x_a ;
  std::vector<double> x_b ;
  std::vector<double> x_g ;

  nbasis = com.nbas() ;
  tmp.resize( 2*nbasis, 2*nbasis) ;
  tmp.setZero() ;
  eulrgrd ( 8, 8, 8, w_a, w_b, w_g, x_a, x_b, x_g, 3) ;
  rotd.fil_mos( nbasis, tmp, 6) ;

  // Loop over the spin integration

  for ( int a=0; a < 8; a++ ){
    for ( int b=0; b < 8; b++ ){
      for ( int g=0; g < 8; g++ ){
        R_s ( com, refd, rotd, x_a[a], x_b[b], x_g[g]) ;
        energy = fockop ( com, H, intarr, refd, rotd, ovl) ;
        w_val = wigner_D ( 0, 0, 0, x_a[a], x_b[b], x_g[g]) ;
        h_k += w_a[a]*w_b[b]*w_g[g]*ovl*energy ;
        n_k += w_a[a]*w_b[b]*w_g[g]*ovl ;
      }
    }
  }
 
  std::cout << "h_k = " << h_k << std::endl ;
  std::cout << "n_k = " << n_k << std::endl ;
  std::cout << "E = " << h_k/n_k << std::endl ;

  return ;
 
} ;

