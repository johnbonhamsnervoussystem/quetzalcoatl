#include<Eigen/Dense>
#include<Eigen/LU>
#include<iostream>
#include "hfwfn.h"
#include "solver.h"
#include "util.h"
/* Utilities that don't belong elsewhere. */

void oao( int nbasis, hfwfn a, Eigen::MatrixXf s) {
  /* Given an overlap and a Slater determinant, check if the determinant is
   * already in the orthogonal atomic orbital basis.  If not, put it into the 
   * oao basis. */
  int iwfnt ;
  float detr ;
  float detru ;
  float d1 = 1e0 ;
  float thresh = 1e-7 ;
  std::complex<float> detc ;
  Eigen::MatrixXf shf_r ; 
  Eigen::MatrixXf mor ; 
  Eigen::MatrixXf tmpr ; 
  Eigen::MatrixXcf shf_c ; 
  Eigen::MatrixXcf moc ; 
  Eigen::MatrixXf tmpc ; 
  
  iwfnt = a.get_wti() ;
  if ( iwfnt % 2 == 1 ){

  /* Real wavefuntion*/

    if ( iwfnt == 1 ){

    /* Restricted wavefuntion*/

      mor.resize( nbasis, nbasis) ;
      tmpr.resize( nbasis, nbasis) ;
      a.get_mos( mor) ;
      tmpr = mor.adjoint()*mor ;
      detr = tmpr.determinant() ;

      if ( abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr = shf_r.adjoint()*s*mor ;
        
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

    } else if ( iwfnt == 3 ){

    /* unrestricted wavefuntion*/

      mor.resize( nbasis, 2*nbasis) ;
      tmpr.resize( nbasis, 2*nbasis) ;
      a.get_mos( mor) ;
      tmpr.block( 0, 0, nbasis, nbasis) = mor.block( 0, 0, nbasis, nbasis).adjoint()*mor.block( 0, 0, nbasis, nbasis) ;
      tmpr.block( 0, nbasis, nbasis, nbasis) = mor.block( 0, nbasis, nbasis, nbasis).adjoint()*mor.block( 0, nbasis, nbasis, nbasis) ;
      detr = tmpr.block( 0, 0, nbasis, nbasis).determinant() ;
      detru = tmpr.block( 0, nbasis, nbasis, nbasis).determinant() ;

      if ( abs(detr - d1) > thresh |  abs(detru - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr.block( 0, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, 0, nbasis, nbasis) ;
        tmpr.block( 0, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

    } else {

    /* Generalized wavefuntion*/

      mor.resize( 2*nbasis, 2*nbasis) ;
      tmpr.resize( 2*nbasis, 2*nbasis) ;
      a.get_mos( mor) ;
      tmpr = mor.adjoint()*mor ;
      detr = tmpr.determinant() ;

      if ( abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr.block( 0, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, 0, nbasis, nbasis) ;
        tmpr.block( 0, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, nbasis, nbasis, nbasis) ;
        tmpr.block( nbasis, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( nbasis, 0, nbasis, nbasis) ;
        tmpr.block( nbasis, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( nbasis, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;


      }

    } else {

  /* complex wavefuntion*/

    if ( iwfnt == 1 ){

    /* Restricted wavefuntion*/

      mor.resize( nbasis, nbasis) ;
      tmpr.resize( nbasis, nbasis) ;
      a.get_mos( mor) ;
      tmpr = mor.adjoint()*mor ;
      detr = tmpr.determinant() ;

      if ( abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr = shf_r.adjoint()*s*mor ;
        
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

    } else if ( iwfnt == 3 ){

    /* unrestricted wavefuntion*/

      mor.resize( nbasis, 2*nbasis) ;
      tmpr.resize( nbasis, 2*nbasis) ;
      a.get_mos( mor) ;
      tmpr.block( 0, 0, nbasis, nbasis) = mor.block( 0, 0, nbasis, nbasis).adjoint()*mor.block( 0, 0, nbasis, nbasis) ;
      tmpr.block( 0, nbasis, nbasis, nbasis) = mor.block( 0, nbasis, nbasis, nbasis).adjoint()*mor.block( 0, nbasis, nbasis, nbasis) ;
      detr = tmpr.block( 0, 0, nbasis, nbasis).determinant() ;
      detru = tmpr.block( 0, nbasis, nbasis, nbasis).determinant() ;

      if ( abs(detr - d1) > thresh |  abs(detru - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr.block( 0, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, 0, nbasis, nbasis) ;
        tmpr.block( 0, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

    } else {

    /* Generalized wavefuntion*/

      mor.resize( 2*nbasis, 2*nbasis) ;
      tmpr.resize( 2*nbasis, 2*nbasis) ;
      a.get_mos( mor) ;
      tmpr = mor.adjoint()*mor ;
      detr = tmpr.determinant() ;

      if ( abs(detr - d1) > thresh ){ 

      /* The determinant is not orthogonal. */

        shf_r.resize( nbasis, nbasis) ;
        shf_c.resize( nbasis, nbasis) ;
        canort( s, shf_c, nbasis) ;
        shf_r = shf_c.real() ;
        tmpr.block( 0, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, 0, nbasis, nbasis) ;
        tmpr.block( 0, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( 0, nbasis, nbasis, nbasis) ;
        tmpr.block( nbasis, 0, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( nbasis, 0, nbasis, nbasis) ;
        tmpr.block( nbasis, nbasis, nbasis, nbasis) = shf_r.adjoint()*s*mor.block( nbasis, nbasis, nbasis, nbasis) ;
  
        a.set_mos( tmpr) ;

      } 

      shf_c.resize( 0, 0) ; 
      shf_r.resize( 0, 0) ;  
      tmpr.resize( 0, 0) ;
      mor.resize( 0, 0) ;

      }  

    }

  }  ;

