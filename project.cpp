#include "constants.h"
#include "common.h"
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <iostream>
#include "evalm.h"
#include "hfwfn.h"
#include "integr.h"
#include "project.h"
#include "tei.h"
#include "util.h"
#include "wigner.h"

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

void phf_drv( common& com, int opt){
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

*/

  if( true ){
     ;
  }

  return ;

} ;

/* 
  Solve for a projected Hartree-Fock wavefunction by repeated 
  diagonalization of the Projected Fock operator.
  
  The theory is outlined in 
  
  Jimenez-hoyos, Carlos A.; Henderson, Thomas M.; Tsuchimochi, T;
  Scuseria, G.; "Projected Hartree-Fock theory" J. Chem. Phys. 
  136, (2012)
*/

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

