#include<Eigen/Dense>
#include<vector>
#include<string>
#include<iostream>
#include "solver.h"
#include "common.h"
#include "hfwfn.h"
#include "hfrout.h"
#include "tei.h"

/* Initialization */
  void hfwfn::init( common& com, std::vector<tei>& intarr, std::string wfn){

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
      wfn_styp = wfn ;
      wfn_ityp = 1 ;

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
      wfn_styp = wfn ;
      wfn_ityp = 2 ;

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
      wfn_styp = wfn ;
      wfn_ityp = 3 ;

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
      wfn_styp = wfn ;
      wfn_ityp = 4 ;

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
      wfn_styp = wfn ;
      wfn_ityp = 5 ;

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
      wfn_styp = wfn ;
      wfn_ityp = 6 ;

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

void hfwfn::fil_mos ( int nbasis, Eigen::Ref<Eigen::MatrixXf> mo, int wfn){

/* Real :: Use the hf class as a container */
  if ( wfn == 1 ) {

    /* Initialize a single block */
    mo_rcof.resize( nbasis, nbasis) ;
    mo_rcof = mo ;

  } else if ( wfn == 3 ) {

    /* Initialize two spin blocks */
    mo_rcof.resize( nbasis, 2*nbasis) ;
    mo_rcof = mo ;

  } else if ( wfn == 5 ) {

    /* Initialize a single block */
    mo_rcof.resize( 2*nbasis, 2*nbasis) ;
    mo_rcof = mo ;

  }

  wfn_ityp = wfn ;

  return ;

} 

void hfwfn::fil_mos ( int nbasis, Eigen::Ref<Eigen::MatrixXcf> mo, int wfn){

/* Complex :: Use the hf class as a container */
  if ( wfn == 2 ) {

    /* Initialize a single block */
    mo_ccof.resize( nbasis, nbasis) ;
    mo_ccof = mo ;

  } else if ( wfn == 4 ) {

    /* Initialize two spin blocks */
    mo_ccof.resize( nbasis, 2*nbasis) ;
    mo_ccof = mo ;

  } else if ( wfn == 6 ) {

    /* Initialize a single block */
    mo_ccof.resize( 2*nbasis, 2*nbasis) ;
    mo_ccof = mo ;

  }

  wfn_ityp = wfn ;

  return ;

} 

/* setter routines */
void hfwfn::set_mos ( Eigen::Ref<Eigen::MatrixXf> mo){
  mo_rcof = mo ;
  return ;
}

void hfwfn::set_mos ( Eigen::Ref<Eigen::MatrixXcf> mo){
  mo_ccof = mo ;
  return ;
}

/* getter routines */
void hfwfn::get_mos ( Eigen::Ref<Eigen::MatrixXf> mo){
  mo = mo_rcof ;
  return ;
}

void hfwfn::get_mos ( Eigen::Ref<Eigen::MatrixXcf> mo){
  mo = mo_ccof ;
  return ;
}

int hfwfn::get_wti ( void){ return wfn_ityp ;}
std::string hfwfn::get_wts ( void){ return wfn_styp ;}


/* Printing routines */
void hfwfn::prt_mos( void) {
    if ( wfn_ityp % 2 == 0 ){
      std::cout << mo_ccof << std::endl ;
    } else {
      std::cout << mo_rcof << std::endl ;
    }
    return ;
}

void hfwfn::prt_eig( void) { std::cout << eig_v << std::endl ; return ; }

void hfwfn::prt_ene( float nn){ std::cout << energy + nn << std::endl ; return ; }

/* Excitation operators */
void hfwfn::ia ( Eigen::Ref<Eigen::MatrixXf> mo, int i, int a){
  /* Return the mo coefficients with orbital i replaced by orbital a */
  mo = mo_rcof ;
  mo.col(i) = mo_rcof.col(a) ;
  mo.col(a) = mo_rcof.col(i) ;
  return ;
}

void hfwfn::ia ( Eigen::Ref<Eigen::MatrixXcf> mo, int i, int a){
  mo = mo_ccof ;
  mo.col(i) = mo_ccof.col(a) ;
  mo.col(a) = mo_ccof.col(i) ;
  return ;
}

