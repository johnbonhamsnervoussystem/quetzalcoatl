#include<Eigen/Dense>
#include<vector>
#include<string>
#include<iostream>
#include "common.h"
#include "hfwfn.h"
#include "hfrout.h"
#include "tei.h"

/* Initialize a hartree-fock solution object whether it is read in or needs to be created. */
  void hfwfn::set( common& com, std::vector<tei>& intarr, std::string wfn){

    int nb ;
    int nele ;
    float e_energy ;
    Eigen::MatrixXf hr ;
    Eigen::MatrixXf sr ;
    Eigen::MatrixXf mor ;
    Eigen::MatrixXcf hc ;
    Eigen::MatrixXcf sc ;
    Eigen::MatrixXcf moc ;
    Eigen::VectorXf eig ;

    if ( wfn == "rrhf" ) {

      /* Initialize a single block */
      mo_rcof.resize( com.nbas(), com.nbas()) ;
      mor.resize( com.nbas(), com.nbas()) ;
      mor.setZero() ;
      eig_v.resize( com.nbas()) ;
      eig.resize( com.nbas()) ;
      hr.resize( com.nbas(), com.nbas()) ;
      sr.resize( com.nbas(), com.nbas()) ;
      hr = com.getH() ;
      sr = com.getS() ;
      wfn_typ = wfn ;

      e_energy = rrhfdia( hr, sr, intarr, com.nbas(), com.nele(), mor, eig) ;

      /* Store the values. */
      mo_rcof = mor ;
      eig_v = eig ;
      energy = e_energy ;

      sr.resize( 0, 0) ;
      hr.resize( 0, 0) ;
      eig.resize( 0) ;
      mor.resize( 0, 0) ;

    } else if ( wfn == "crhf" ) {

      /* Initialize a single block */
      mo_ccof.resize( com.nbas(), com.nbas()) ;
      moc.resize( com.nbas(), com.nbas()) ;
      moc.setZero() ;
      eig_v.resize( com.nbas()) ;
      eig.resize( com.nbas()) ;
      hc.resize( com.nbas(), com.nbas()) ;
      sc.resize( com.nbas(), com.nbas()) ;
      hc.real() = com.getH() ;
      sc.real() = com.getS() ;
      wfn_typ = wfn ;

      e_energy = crhfdia( hc, sc, intarr, com.nbas(), com.nele(), moc, eig) ;

      /* Store the values. */
      mo_ccof = moc ;
      eig_v = eig ;
      energy = e_energy ;

      sc.resize( 0, 0) ;
      hc.resize( 0, 0) ;
      eig.resize( 0) ;
      mor.resize( 0, 0) ;

    } else if ( wfn == "ruhf" ) {

      /* Initialize a single block */
      mo_rcof.resize( com.nbas(), 2*com.nbas()) ;
      mor.resize( com.nbas(), 2*com.nbas()) ;
      mor.setZero() ;
      eig_v.resize( 2*com.nbas()) ;
      eig.resize( 2*com.nbas()) ;
      hr.resize( com.nbas(), com.nbas()) ;
      sr.resize( com.nbas(), com.nbas()) ;
      hr = com.getH() ;
      sr = com.getS() ;
      wfn_typ = wfn ;

      e_energy = ruhfdia( hr, sr, intarr, com.nbas(), com.nalp(), com.nbet(), mor.block( 0, 0, com.nbas(), com.nbas()), mor.block( 0, com.nbas(), com.nbas(), com.nbas()), eig) ;

      /* Store the values. */
      mo_rcof = mor ;
      eig_v = eig ;
      energy = e_energy ;

      sr.resize( 0, 0) ;
      hr.resize( 0, 0) ;
      eig.resize( 0) ;
      mor.resize( 0, 0) ;

    } else if ( wfn == "cuhf" ) {

      /* Initialize a single block */
      mo_ccof.resize( com.nbas(), 2*com.nbas()) ;
      moc.resize( com.nbas(), 2*com.nbas()) ;
      moc.setZero() ;
      eig_v.resize( 2*com.nbas()) ;
      eig.resize( 2*com.nbas()) ;
      hc.resize( com.nbas(), com.nbas()) ;
      sc.resize( com.nbas(), com.nbas()) ;
      hc.setZero() ;
      sc.setZero() ;
      hc.real() = com.getH() ;
      sc.real() = com.getS() ;
      wfn_typ = wfn ;

      std::cout << "Just before complex unrestricted " << std::endl ;
      e_energy = cuhfdia( hc, sc, intarr, com.nbas(), com.nalp(), com.nbet(), moc.block( 0, 0, com.nbas(), com.nbas()), moc.block( 0, com.nbas(), com.nbas(), com.nbas()), eig) ;

      /* Store the values. */
      mo_ccof = moc ;
      eig_v = eig ;
      energy = e_energy ;

      sc.resize( 0, 0) ;
      hc.resize( 0, 0) ;
      eig.resize( 0) ;
      moc.resize( 0, 0) ;

    } else if ( wfn == "rghf" ) {

      /* Initialize a single block */
      mo_rcof.resize( 2*com.nbas(), 2*com.nbas()) ;
      mor.resize( 2*com.nbas(), 2*com.nbas()) ;
      mor.setZero() ;
      eig_v.resize( 2*com.nbas()) ;
      eig.resize( 2*com.nbas()) ;
      hr.resize( com.nbas(), com.nbas()) ;
      sr.resize( com.nbas(), com.nbas()) ;
      hr = com.getH() ;
      sr = com.getS() ;
      wfn_typ = wfn ;

      e_energy = rghfdia( hr, sr, intarr, com.nbas(), com.nele(), mor, eig) ;

      /* Store the values. */
      mo_rcof = mor ;
      eig_v = eig ;
      energy = e_energy ;

      sr.resize( 0, 0) ;
      hr.resize( 0, 0) ;
      eig.resize( 0) ;
      mor.resize( 0, 0) ;

    } else if ( wfn == "cghf" ) {

      /* Initialize a single block */
      mo_ccof.resize( 2*com.nbas(), 2*com.nbas()) ;
      moc.resize( 2*com.nbas(), 2*com.nbas()) ;
      moc.setZero() ;
      eig_v.resize( 2*com.nbas()) ;
      eig.resize( 2*com.nbas()) ;
      hc.resize( com.nbas(), com.nbas()) ;
      sc.resize( com.nbas(), com.nbas()) ;
      hc.setZero() ;
      sc.setZero() ;
      hc.real() = com.getH() ;
      sc.real() = com.getS() ;
      wfn_typ = wfn ;

      e_energy = cghfdia( hc, sc, intarr, com.nbas(), com.nele(), moc, eig) ;

      /* Store the values. */
      mo_ccof = moc ;
      eig_v = eig ;
      energy = e_energy ;

      sc.resize( 0, 0) ;
      hc.resize( 0, 0) ;
      eig.resize( 0) ;
      moc.resize( 0, 0) ;

    }

    return ;

} ;

/* Printing routines */

void hfwfn::prt_mos( void) {
    if ( wfn_typ[0] == 'c' ){
      std::cout << mo_ccof << std::endl ;
    } else {
      std::cout << mo_rcof << std::endl ;
    }
    return ;
}

void hfwfn::prt_eig( void) { std::cout << eig_v << std::endl ; return ; }

void hfwfn::prt_ene( float nn){ std::cout << energy + nn << std::endl ; return ; }

//  void hfwfn::set( common& com, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> mocof, Eigen::Ref<Eigen::VectorXf> eig, float e_energy, std::string wfn){
//    bool init=false ;
//    Eigen::MatrixXf h ;
//    Eigen::MatrixXf s ;
//
//    init = eig.isZero(0) ;
//
//    if ( wfn == "rhf" ){
//
//      /* Initialize a single block */
//      mo_rcof.resize( com.nbas(), com.nbas()) ;
//      eig_v.resize( com.nbas()) ;
//      eig.resize( com.nbas()) ;
//      wfn_typ = "rrhf" ;
//
//      if ( init ){
//
//        /* The vecctor of eigenvalues was empty so we solve for a wavefunction here. */
//        h.resize( com.nbas(), com.nbas()) ;
//        s.resize( com.nbas(), com.nbas()) ;
//        h = com.getH() ;
//        s = com.getS() ;
//
//        if ( zero_den ) {
//
//          /* No density was read in so we don't have an initial guess. */
//          mocof.resize( com.nbas(), com.nbas()) ;
//          mocof.setZero() ;
//
//        }
//       
//        /* Find a self-consistent wavefunction */
//        e_energy = rrhfdia( h, s, intarr, com.nbas(), com.nele(), mocof, eig) ;
//
//      }
//
//      /* Store the values. */
//      mo_rcof = mocof ;
//      eig_v = eig ;
//      energy = e_energy ;
//
//    } else if ( wfn == "uhf" ){
//
//      /* Initialize a double long block */
//      mo_rcof.resize( com.nbas(), 2*com.nbas()) ;
//      eig.resize( 2*com.nbas()) ;
//      wfn_typ = "ruhf" ;
//
//      if ( init ){
//
//        /* The vecctor of eigenvalues was empty so we solve for a wavefunction here. */
//        h.resize( com.nbas(), com.nbas()) ;
//        s.resize( com.nbas(), com.nbas()) ;
//        h = com.getH() ;
//        s = com.getS() ;
//
//        if ( zero_den ) {
//
//          /* No density was read in so we don't have an initial guess. */
//          mocof.resize( com.nbas(), 2*com.nbas()) ;
//          mocof.setZero() ;
//
//        }
//       
//        /* Find a self-consistent wavefunction */
//        e_energy = ruhfdia( h, s, intarr, com.nbas(), com.nalp(), com.nbet(), mocof.block( 0, 0, com.nbas(), com.nbas()), mocof.block( 0, com.nbas(), com.nbas(), com.nbas()), eig) ;
//
//      }
//
//      /* Store the values. */
//      mo_rcof = mocof ;
//      eig_v = eig ;
//      energy = e_energy ;
//
//    } else if ( wfn == "ghf" ){
//
//      /* Initialize a single block */
//      mo_rcof.resize( 2*com.nbas(), 2*com.nbas()) ;
//      eig.resize( 2*com.nbas()) ;
//      wfn_typ = "rghf" ;
//
//      if ( init ){
//
//        /* The vecctor of eigenvalues was empty so we solve for a wavefunction here. */
//        h.resize( com.nbas(), com.nbas()) ;
//        s.resize( com.nbas(), com.nbas()) ;
//        h = com.getH() ;
//        s = com.getS() ;
//
//        if ( zero_den ) {
//
//          /* No density was read in so we don't have an initial guess. */
//          mocof.resize( 2*com.nbas(), 2*com.nbas()) ;
//          mocof.setZero() ;
//
//        }
//       
//        /* Find a self-consistent wavefunction */
//        e_energy = rghfdia( h, s, intarr, com.nbas(), com.nele(), mocof, eig) ;
//
//      }
//
//      /* Store the values. */
//      mo_rcof = mocof ;
//      eig_v = eig ;
//      energy = e_energy ;
//
//    } 
//
//    return ;
//
//} ;
//
//  void hfwfn::set( common& com, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> mocof, Eigen::Ref<Eigen::VectorXf> eig, float e_energy, std::string wfn){
//    bool init=false ;
//    bool zero_den=true ;
//    Eigen::MatrixXcf h ;
//    Eigen::MatrixXcf s ;
//
//    init = eig.isZero(0) ;
//    zero_den = mocof.isZero(0) ;
//
//    if ( wfn == "rhf" ){
//
//      /* Initialize a single block */
//      mo_ccof.resize( com.nbas(), com.nbas()) ;
//      eig.resize( com.nbas()) ;
//      wfn_typ = "crhf" ;
//
//      if ( init ){
//
//        /* The vecctor of eigenvalues was empty so we solve for a wavefunction here. */
//        h.resize( com.nbas(), com.nbas()) ;
//        s.resize( com.nbas(), com.nbas()) ;
//        h.real() = com.getH() ;
//        s.real() = com.getS() ;
//
//        if ( zero_den ) {
//
//          /* No density was read in so we don't have an initial guess. */
//          mocof.resize( com.nbas(), com.nbas()) ;
//          mocof.setZero() ;
//
//        }
//       
//        /* Find a self-consistent wavefunction */
//        e_energy = crhfdia( h, s, intarr, com.nbas(), com.nele(), mocof, eig) ;
//
//      }
//
//      /* Store the values. */
//      mo_ccof = mocof ;
//      eig_v = eig ;
//      energy = e_energy ;
//
//    } else if ( wfn == "uhf" ){
//
//      /* Initialize a double long block */
//      mo_ccof.resize( com.nbas(), 2*com.nbas()) ;
//      eig.resize( 2*com.nbas()) ;
//      wfn_typ = "cuhf" ;
//
//      if ( init ){
//
//        /* The vecctor of eigenvalues was empty so we solve for a wavefunction here. */
//        h.resize( com.nbas(), com.nbas()) ;
//        s.resize( com.nbas(), com.nbas()) ;
//        h.real() = com.getH() ;
//        s.real() = com.getS() ;
//
//        if ( zero_den ) {
//
//          /* No density was read in so we don't have an initial guess. */
//          mocof.resize( com.nbas(), 2*com.nbas()) ;
//          mocof.setZero() ;
//
//        }
//       
//        /* Find a self-consistent wavefunction */
//        e_energy = cuhfdia( h, s, intarr, com.nbas(), com.nalp(), com.nbet(), mocof.block( 0, 0, com.nbas(), com.nbas()), mocof.block( 0, com.nbas(), com.nbas(), com.nbas()), eig) ;
//
//      }
//
//      /* Store the values. */
//      mo_ccof = mocof ;
//      eig_v = eig ;
//      energy = e_energy ;
//
//    } else if ( wfn == "ghf" ){
//
//      /* Initialize a single block */
//      mo_ccof.resize( 2*com.nbas(), 2*com.nbas()) ;
//      eig.resize( 2*com.nbas()) ;
//      wfn_typ = "cghf" ;
//
//      if ( init ){
//
//        /* The vecctor of eigenvalues was empty so we solve for a wavefunction here. */
//        h.resize( com.nbas(), com.nbas()) ;
//        s.resize( com.nbas(), com.nbas()) ;
//        h.real() = com.getH() ;
//        s.real() = com.getS() ;
//
//        if ( zero_den ) {
//
//          /* No density was read in so we don't have an initial guess. */
//          mocof.resize( 2*com.nbas(), 2*com.nbas()) ;
//          mocof.setZero() ;
//
//        }
//       
//        /* Find a self-consistent wavefunction */
//        e_energy = cghfdia( h, s, intarr, com.nbas(), com.nele(), mocof, eig) ;
//
//      }
//
//      /* Store the values. */
//      mo_ccof = mocof ;
//      eig_v = eig ;
//      energy = e_energy ;
//
//    } 
//
//    return ;
//
//} ;
//
