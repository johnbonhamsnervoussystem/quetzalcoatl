#include<cmath>
#include<iostream>
#include<complex>
#include<vector>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include "solver.h"
#include "evalm.h"

void canort( Eigen::Ref<Eigen::MatrixXf> const s, Eigen::Ref<Eigen::MatrixXcf> xs , int& dim) {
  /*  
 *    Find an orthogonal transformation using canonical orthogonalization for a real matrix
 *    */

  int indx=-1 ;
  float thresh = 1e-7 ;
  std::complex<float> eig ;
  Eigen::EigenSolver<Eigen::MatrixXf> s_diag( s, true) ;
  Eigen::MatrixXcf s_U ;
  Eigen::VectorXcf s_eig ;

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

void canort( Eigen::Ref<Eigen::MatrixXcf> const s, Eigen::Ref<Eigen::MatrixXcf> xs , int& dim) {
  /*  
 *    Find an orthogonal transformation using canonical orthogonalization for a complex matrix
 *    */

  int indx=-1 ;
  float thresh = 1e-7 ;
  std::complex<float> eig ;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcf> s_diag( s, true) ;
  Eigen::MatrixXcf s_U ;
  Eigen::VectorXcf s_eig ;

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

