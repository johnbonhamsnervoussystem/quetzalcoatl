#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "numer.h"
#include <vector>

void numerical_derivative( double& h, Eigen::Ref<Eigen::MatrixXd> S, double (*f)( int&, cd&, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>), int& nbas, cd& nrm, Eigen::Ref<Eigen::MatrixXd> T, Eigen::Ref<Eigen::MatrixXd> U) {
/*
  I can generalize this routine later.  For now let's use it for what is needed.
  h - the step size to be used
  S - the matrix to differentiate w.r.t.
  (*f) - the function to evaluate 
  T - additional arguments to (*f)
*/
  double t, w ;
  std::vector<double> v ;
  Eigen::MatrixXd::Index r = S.rows() ;
  Eigen::MatrixXd::Index c = S.cols() ;
  Eigen::MatrixXd Sx ( r, c) ;

  for( int i = 0; i < r; i++){
    for( int j = 0; j < c; j++){
      t = S( i, j) ;
      v.clear() ;
      for( int k = 0; k < 4; k++){
        switch( k){
          case 0 : S( i, j) = t + h + h ;
            break ;
          case 1 : S( i, j) = t + h ;
            break ;
          case 2 : S( i, j) = t - h ;
            break ;
          case 3 : S( i, j) = t - h - h ;
            break ;
          }
        w = (*f)( nbas, nrm, S, T, U) ;
        v.push_back( w) ;
        }
      Sx( i, j) = five_point( h, v) ;
      S( i, j) = t ;
      }
    }

  std::cout << Sx << std::endl ;

  return ;

} ;

void numerical_derivative( double& h, Eigen::Ref<Eigen::MatrixXcd> S, cd (*f)( int&, cd&, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>), int& nbas, cd& nrm, Eigen::Ref<Eigen::MatrixXcd> T, Eigen::Ref<Eigen::MatrixXcd> U) {
/*
  I can generalize this routine later.  For now let's use it for what is needed.
  h - the step size to be used
  S - the matrix to differentiate w.r.t.
  (*f) - the function to evaluate 
  T - additional arguments to (*f)

  Complex version of the numerical derivatives
  Do the real part and then the complex part.
*/

  double t, t2 ;
  cd w ;
  std::vector<cd> v ;
  Eigen::MatrixXcd::Index r = S.rows() ;
  Eigen::MatrixXcd::Index c = S.cols() ;
  Eigen::MatrixXcd Sx ( r, c) ;
  Sx.setZero() ;

/*
  Real
*/

  std::cout << " h in numer " << h << std::endl ;

  for( int i = 0; i < r; i++){
    for( int j = i; j < c; j++){
      t = S( i, j).real() ;
      t2 = S( i, j).imag() ;
      v.clear() ;
      for( int k = 0; k < 4; k++){
        switch( k){
          case 0 : S( i, j) = cd( t + h + h, t2) ;
            break ;
          case 1 : S( i, j) = cd( t + h, t2) ;
            break ;
          case 2 : S( i, j) = cd( t - h, t2) ;
            break ;
          case 3 : S( i, j) = cd( t - h - h, t2) ;
            break ;
          }
        w = (*f)( nbas, nrm, S, T, U) ;
        v.push_back( w) ;
        }
      Sx( i, j) += five_point( h, v)/z2 ;
      S( i, j) = cd( t, t2) ;
      }
    }

/*
  Complex
*/

  for( int i = 0; i < r; i++){
    for( int j = i; j < c; j++){
      t = S( i, j).imag() ;
      t2 = S( i, j).real() ;
      v.clear() ;
      for( int k = 0; k < 4; k++){
        switch( k){
          case 0 : S( i, j) = cd( t2, t + h + h) ;
            break ;
          case 1 : S( i, j) = cd( t2, t + h) ;
            break ;
          case 2 : S( i, j) = cd( t2, t - h) ;
            break ;
          case 3 : S( i, j) = cd( t2, t - h - h) ;
            break ;
          }
        w = (*f)( nbas, nrm, S, T, U) ;
        v.push_back( w) ;
        }
      Sx( i, j) -= zi*five_point( h, v)/z2 ;
      S( i, j) = cd( t2, t) ;
      }
    }

  std::cout << Sx << std::endl ;

  Sx.resize( 0, 0) ;

  return ;

} ;

double five_point( const double& h, std::vector<double>& f){
  double d12 = 12.0e0 ;

  if ( f.size() != 4 ){
    std::cout << " Not enough points passed to five point " << std::endl ;
    throw std::exception() ;
    } 

  return (-f[0] + d8*f[1] - d8*f[2] + f[3])/( d12*h) ;

} ;

cd five_point( const double& h, std::vector<cd>& f){
  cd dh ( h, d0) ;
  cd z8 ( d8, d0) ;
  cd z12 ( 12.0e0, d0) ;

  if ( f.size() != 4 ){
    std::cout << " Not enough points passed to five point " << std::endl ;
    throw std::exception() ;
    } 

  return (-f[0] + z8*f[1] - z8*f[2] + f[3])/( z12*dh) ;

} ;

