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


double rhfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nele, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, bool lshift, double lvlshft, const int& maxit, const double& thresh) {

  /* 
    Restricted Hartree-Fock solved by repeated diagonalization.
  */
  matrix f, g, p, p_prev ;
  Eigen::SelfAdjointEigenSolver<matrix> f_diag ;
  int iter=0 ;
  int occ ;
  typename matrix::Scalar energy, t, shift, two ;
  time_dbg rhfdia_time = time_dbg("rhfdia") ;

  two = static_cast<typename matrix::Scalar>(d2) ;
  occ = nele/2 ;
  if (lshift){
    shift = static_cast<typename matrix::Scalar>(lvlshft) ;
    }
  f.resize( nbasis, nbasis) ;
  g.resize( nbasis, nbasis) ;
  p.resize( nbasis, nbasis) ;
  p_prev.resize( nbasis, nbasis) ;

  /* If c has something in it use it as the initial guess. */
  if( c.isZero(0) ) {
    f = h ;
  } else {
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    W->contract(p) ;
    g = W->getG() ;
    f = h + g ;
  }

  while ( iter++ < maxit ) {
    p_prev = p ;
    if ( lshift && iter > 3){
      f += -shift*p ;
      }
    f_diag.compute(f) ;
    c = f_diag.eigenvectors() ;
    p = c.block( 0, 0, nbasis, occ)*c.block( 0, 0, nbasis, occ).adjoint() ;
    W->contract( p) ;
    g = W->getG() ;
    f = h + g ;
    g = p*( h + f) ;

    energy = g.trace() ;

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
