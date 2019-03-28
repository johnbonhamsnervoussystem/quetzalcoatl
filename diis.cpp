#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include "constants.h"
#include "diis.h"
#include "qtzcntrl.h"
#include "qtzio.h"

template< class s>
void diis<s>::update( Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> f, bool& ex_cond){

  int i, j, np ;
  /*
    Calculate the error vector fp - pf and save in e_diis
  */

  scr = f*p - p*f ;
  t = Eigen::Map<Eigen::MatrixXd>( scr.data(), scr.cols()*scr.rows(), 1) ;
  e_diis.push_back( t) ;

  /*
    Store the current fock matrix
  */
  t = Eigen::Map<Eigen::MatrixXd>( f.data(), f.cols()*f.rows(), 1) ;
  f_diis.push_back( t) ;

  np = e_diis.size() ;

  /*
    Check if we should do the DIIS
  */
  if ( ! diis_on){
    if ( scr.cwiseAbs().maxCoeff() < 0.1 && ex_cond ){
      std::cout << " Maximum error is below threshhold. Turning on diis" << std::endl ;
      diis_on = true ;
      }
    }

  B.setZero() ;

  if ( np > n_e && diis_on){
    B.block(0, 0, np, np).setConstant( d1) ;
    B( np-1, np-1) = d0 ;
    for ( i = 0; i < np - 1; i++){
      t = e_diis[i+1] - e_diis[i] ;
      for ( j = i; j < np - 1; j++){
        t1 = e_diis[j+1] - e_diis[j] ;
        B( i, j) = t.dot( t1) ;
        B( j, i) = B( i, j) ;
        }
      }
    v.setZero() ;
    c.setZero() ;
    v( np-1) = d1 ;
    c = B.colPivHouseholderQr().solve(v) ;

    t.setZero() ;

    for( i = 0; i < np-1; i++){
      t += c(i)*f_diis[i] ;
      }

    /*
      t is our new density.  Replace the old guess and send it back into the world.
    */
    f_diis.back() = t ;
    f = Eigen::Map<Eigen::MatrixXd>( t.data(), f.rows(), f.cols()) ;
    }

  while ( np > n_e ){
    e_diis.erase( e_diis.begin()) ;
    f_diis.erase( f_diis.begin()) ;
    np = e_diis.size() ;
    }

  return ;

} ;

template void diis<double>::update( Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> f, bool& ex_cond) ;

template< class s>
void diis<s>::update( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> f, bool& ex_cond){

  int i, j, np ;
  /*
    Calculate the error vector fp - pf and save in e_diis
  */

  scr = f*p - p*f ;
  t = Eigen::Map<Eigen::MatrixXcd>( scr.data(), scr.cols()*scr.rows(), 1) ;
  e_diis.push_back( t) ;

  /*
    Store the current fock matrix
  */
  t = Eigen::Map<Eigen::MatrixXcd>( f.data(), f.cols()*f.rows(), 1) ;
  f_diis.push_back( t) ;

  np = e_diis.size() ;

  /*
    Check if we should do the DIIS
  */

  if ( ! diis_on){
    if ( scr.cwiseAbs().maxCoeff() < 0.1 && ex_cond ){
      std::cout << " Maximum error is below threshhold. Turning on diis" << std::endl ;
      diis_on = true ;
      }
    }

  B.setZero() ;

  if ( np > n_e && diis_on){
    B.row( 2*n_e).setConstant( -z1) ;
    B.col( 2*n_e).setConstant( -z1) ;
    B.row( 2*n_e+1).head(n_e).setConstant( -zi) ;
    B.row( 2*n_e+1).tail(n_e+2).setConstant( zi) ;
    B.col( 2*n_e+1).head(n_e).setConstant( zi) ;
    B.col( 2*n_e+1).tail(n_e+2).setConstant( -zi) ;
    B.block( 2*n_e, 2*n_e, 2, 2).setConstant( z0) ;
    for ( i = 0; i < np-1; i++){
      t = e_diis[i+1].conjugate() - e_diis[i].conjugate() ;
      for ( j = i; j < np - 1; j++){
        t1 = e_diis[j+1] - e_diis[j] ;
        B( i, j) = t.dot( t1) ;
        B( j, i) = std::conj(B( i, j)) ;
        }
      }
    B.block( n_e, n_e, n_e, n_e) = B.block( 0, 0, n_e, n_e).conjugate() ;
    v.setZero() ;
    c.setZero() ;
    v( 2*n_e) = -z2 ;
    c = B.colPivHouseholderQr().solve(v) ;

    t.setZero() ;

    for( i = 0; i < n_e; i++){
      t += c(i)*f_diis[i] ;
      }

    /*
      t is our new fock matrix.  Replace the old fock and send the new one out into the world.
    */
    f_diis.back() = t ;
    f = Eigen::Map<Eigen::Matrix< s, Eigen::Dynamic, Eigen::Dynamic>>( t.data(), f.rows(), f.cols()) ;
    }

  while ( np > n_e ){
    e_diis.erase( e_diis.begin()) ;
    f_diis.erase( f_diis.begin()) ;
    np = e_diis.size() ;
    }

  return ;

} ;

template void diis<std::complex<double>>::update( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> f, bool& ex_cond) ;

void diis_control::set_switch( int opt){

  control_switch = opt ;

  return ;

} ;

void diis_control::toggle( double e){

  /* Check if diis should be turned on */
  if ( ! diis_switch){
    if ( control_switch == 1) {
      /* Check against the previous value and then update the previous value */
      if ( std::abs(e - cntl_ref) < ediff_thr){
        diis_switch = true ;
        }
      cntl_ref = e ;
    } else if ( control_switch == 2) {
      /* Check that it is below a certain value by the threshhold */
      if ( e - cntl_ref < -ediff_thr){
        diis_switch = true ;
        }
    } else {
      qtzcntrl::shutdown( "DIIS is requrested but no condition set to turn it on.") ;
      }
    }

  return ;

} ;
