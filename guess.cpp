#include "constants.h"
#include <iostream>
#include <Eigen/Dense>
#include "guess.h"
#include "qtzio.h"
#include "wfn.h"

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k) {
/*
  Given an initial real Slater Determinant from Hartree-Fock, return
  a thermalized rho and kappa guess.
*/
  double Ef, djunk ;
  Eigen::MatrixXd A, B, C ;
  wfn < double, Eigen::Dynamic, Eigen::Dynamic> w ;
  A.resize( nbas, nbas) ;
  B.resize( nbas, nbas) ;
  C.resize( nbas, nbas) ;
  A.setZero() ;
  B.setZero() ;
  C.setZero() ;
  w.moc.resize( nbas, nbas) ;
  w.eig.resize( nbas) ;
  load_wfn( w) ;
  Ef = (w.eig( nele - 1) + w.eig( nele))/d2 ;
  for( int i=0; i < nbas; i++){
    djunk = (w.eig( i) - Ef)/(kb*3.0e4) ;
    A( i, i) = d1/( d1 + std::exp( djunk)) ;
    B( i, i) = std::sqrt( A( i, i) - A( i, i)*A(i, i)) ;
    }

  C = w.moc*A*w.moc.adjoint() ;
  print_mat( C, " initial rho") ;

  if ( p.rows() == nbas ){
    p.real() = C ;
  } else {
    p.block( 0, 0, nbas, nbas).real() = C ;
    p.block( nbas, nbas, nbas, nbas).real() = C ;
    }

  C = w.moc*B*w.moc.adjoint() ;
  print_mat( C, " initial kappa") ;

  if ( k.rows() == nbas ){
    k.real() = C ;
  } else {
    k.block( 0, nbas, nbas, nbas).real() = C ;
    k.block( nbas, 0, nbas, nbas).real() = -C ;
    }

  w.eig.resize( 0) ;
  w.moc.resize( 0, 0) ;
  C.resize( 0, 0) ;
  B.resize( 0, 0) ;
  A.resize( 0, 0) ;

  return ;

}

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> k) {
/*
  Given an initial real Slater Determinant from Hartree-Fock, return
  a thermalized rho and kappa guess.
*/
  double Ef, djunk ;
  Eigen::MatrixXd A, B, C ;
  wfn < double, Eigen::Dynamic, Eigen::Dynamic> w ;
  A.resize( nbas, nbas) ;
  B.resize( nbas, nbas) ;
  C.resize( nbas, nbas) ;
  A.setZero() ;
  B.setZero() ;
  C.setZero() ;
  w.moc.resize( nbas, nbas) ;
  w.eig.resize( nbas) ;
  load_wfn( w) ;
  Ef = (w.eig( nele - 1) + w.eig( nele))/d2 ;
  for( int i=0; i < nbas; i++){
    djunk = (w.eig( i) - Ef)/(kb*3.0e4) ;
    A( i, i) = d1/( d1 + std::exp( djunk)) ;
    B( i, i) = std::sqrt( A( i, i) - A( i, i)*A(i, i)) ;
    }

  C = w.moc*A*w.moc.adjoint() ;
  print_mat( C, " initial rho") ;

  if ( p.rows() == nbas ){
    p = C ;
  } else {
    p.block( 0, 0, nbas, nbas) = C ;
    p.block( nbas, nbas, nbas, nbas) = C ;
    }

  C = w.moc*B*w.moc.adjoint() ;
  print_mat( C, " initial kappa") ;

  if ( k.rows() == nbas ){
    k = C ;
  } else {
    k.block( 0, nbas, nbas, nbas) = C ;
    k.block( nbas, 0, nbas, nbas) = -C ;
    }

  w.eig.resize( 0) ;
  w.moc.resize( 0, 0) ;
  C.resize( 0, 0) ;
  B.resize( 0, 0) ;
  A.resize( 0, 0) ;

  return ;

}

