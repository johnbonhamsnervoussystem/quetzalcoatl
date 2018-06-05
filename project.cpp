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
 *
 * Wish list -
 *   -HF
 *     - Spin
 *     - C.C.
 *     - Point Group
 *   -HFB
 *     - All of the above
 *     - Number
 *
 * */

/* 
 * Solve for a projected Hartree-Fock wavefunction by repeated 
 * diagonalization of the Projected Fock operator.
 *
 * The theory is outlined in 
 *
 * Jimenez-hoyos, Carlos A.; Henderson, Thomas M.; Tsuchimochi, T;
 * Scuseria, G.; "Projected Hartree-Fock theory" J. Chem. Phys. 
 * 136, (2012)
 *
 * */

void e_phf( common& com, hfwfn& refd, Eigen::Ref<Eigen::MatrixXcf> H, std::vector<tei>& intarr) {

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
  float Intg = 0.0 ;
  cf w_val ;
  cf energy ;
  cf ovl ;
  cf n_k = cf( 0.0, 0.0) ;
  cf h_k = cf( 0.0, 0.0) ;
  Eigen::MatrixXcf tmp ;

 // One reference determinant. One for storage space.

  hfwfn rotd ;
 
 // integration grid.

  std::vector<float> w_a ;
  std::vector<float> w_b ;
  std::vector<float> w_g ;
  std::vector<float> x_a ;
  std::vector<float> x_b ;
  std::vector<float> x_g ;

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

