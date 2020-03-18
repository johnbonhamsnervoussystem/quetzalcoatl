#include <algorithm>
#include <cmath>
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include "hfwfn.h"
#include "integr.h"
#include <iostream>
#include "solver.h"
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include <vector>
#include "wigner.h"

/* Utilities that don't belong elsewhere. */

  void bonnet_r ( int n, double x, double& p2) {
  /* 
    Use Bonnet's recursion formula to evaluate the nth Legendre Polynomial
    at the point x.
   
    P_{m+1}(x) = ((2n + 1)xP_{m}(x) - nP_{n-1}(x))/(m + 1)
  */
    double p0 ;
    double p1 = 1.0 ;
  
    // Handle edge cases and
    if ( n < 0 ){
     std::cout << "Warning: Negative polynomial in bonnet_r" << std::endl ;
     return ;
    }
  
    if ( n == 0 ){ p2 = 1.0; return ; }
    if ( n == 1 ){ p2 = x ; return ; }
  
    p2 = x ;
   
    for ( int i = 1; i < n ; i++){
      p0 = p1 ;
      p1 = p2 ;
      p2 = (static_cast<double>(2*i + 1)*x*p1 - i*p0)/static_cast<double>(i + 1) ;
    }
   
    return ;
  
  }
  
  void bonnet_r( int n, double x, double& p2, double& dp) {
  /* 
    Use Bonnet's recursion formula to evaluate the nth Legendre Polynomial
    at the point x.  Additionally, return the derivative of P_{n}
   
    P_{m+1}(x) = ((2n + 1)xP_{m}(x) - nP_{n-1}(x))/(m + 1)
  */ 
    double p0 ;
    double p1 = 1.0 ;
  
    // Handle edge cases and
    if ( n < 0 ){
     std::cout << "Warning: Negative polynomial in bonnet_r" << std::endl ;
     return ;
    }
  
    if ( n == 0 ){ p2 = 1.0; dp=d0; return ; }
    if ( n == 1 ){ p2 = x; dp=1.0; return ; }
  
    p2 = x ;
   
    for ( int i = 1; i < n ; i++){
      p0 = p1 ;
      p1 = p2 ;
      p2 = (static_cast<double>(2*i + 1)*x*p1 - i*p0)/static_cast<double>(i + 1) ;
    }
  
    dp = static_cast<double>(n)*(x*p2 - p1)/(x*x - 1.0) ;
  
    return ;
  
  }


  double fact ( int n) {
    /* Return the factorial of n */
    double f = d1 ;
  
    if ( n <= 1 ) { return f ;}
  
    while ( n > 1 ) {
      f = f*static_cast<double>(n) ;
      n+= -1 ;
    }
  
    return f ;
  
  } ;

  double factfact ( int n) {
    /* Return the double factorial of n */
    double f=1.0 ;
  
    if ( n < 1 ) { return 1.0 ;}
  
    while ( n > 1 ) {
      f = f*double(n) ;
      n+= -2 ;
    }
  
    return f ;
  
  } ;

  double fboys( int k, double t) {
      double f, m, n ;
  /*
  Let's not mess around.  Just compute this
  */
  
      m = static_cast<double>(k) + d1/d2 ;
      n = static_cast<double>(k) + 3.0e0/d2 ;
      f = gsl_sf_hyperg_1F1( m, n, -t)/( d2*static_cast<double>(k) + d1) ;
      
      return f ;
  
    } ;
