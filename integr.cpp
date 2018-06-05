#include "constants.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "integr.h"
/* Routines to set up numerical integration grids */

void bonnet_r ( int n, float x, float& p2) {
/* Use Bonnet's recursion formula to evaluate the nth Legendre Polynomial
 * at the point x.
 *
 * P_{m+1}(x) = ((2n + 1)xP_{m}(x) - nP_{n-1}(x))/(m + 1)
 * */
  float p0 ;
  float p1 = 1.0 ;

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
    p2 = (static_cast<float>(2*i + 1)*x*p1 - i*p0)/static_cast<float>(i + 1) ;
  }
 
  return ;

}

void bonnet_r( int n, float x, float& p2, float& dp) {
/* Use Bonnet's recursion formula to evaluate the nth Legendre Polynomial
 * at the point x.  Additionally, return the derivative of P_{n}
 *
 * P_{m+1}(x) = ((2n + 1)xP_{m}(x) - nP_{n-1}(x))/(m + 1)
 * */
  float p0 ;
  float p1 = 1.0 ;

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
    p2 = (static_cast<float>(2*i + 1)*x*p1 - i*p0)/static_cast<float>(i + 1) ;
  }

  dp = static_cast<float>(n)*(x*p2 - p1)/(x*x - 1.0) ;

  return ;

}


void gauleg ( float lv, float uv, int n, std::vector<float>& x, std::vector<float>& w) {
/*
 * This is a routine to set up a Gauss-Legendre quadrature grid.
 * This code is adapted from Carlos Jimenez-Hoyos who obtained
 * the original implementation from
 *
 * Teukolsky, Vetterling and Flannery, Numerical Recipes in F77
 * 2nd Ed. 2005
 *
 * lv - lower value of the integration grid
 * uv - upper value of the integration grid
 * n - number of grid points
 * x - abscissas of the quadrature
 * w - weights of the quadrature
 *
 * */

  int m ;
  float z ;
  float z1 ;
  float xm ;
  float xl ;
  float py ;
  float dpy ;
  float eps = 1.0e-7 ;

  m = (n + 1)/2 ;
  xm = (uv + lv)/2.0 ;
  xl = (uv - lv)/2.0 ;

  for ( int i=1; i <=m; i++ ) {

    z = cos ( pi*(static_cast<float>(i)-0.25)/(static_cast<float>(n)+0.5)) ;

    z1 = z + 1.0 ;
    dpy = 1.0 ;

    while ( std::abs(z-z1) > eps ) {

      bonnet_r ( n, z, py, dpy) ;
      z1 = z ;

      z = z1 - py/dpy ;
    }

    x[i-1] = xm - xl*z ;
    x[n-i] = xm + xl*z ;

    w[i-1] = 2.0*xl/((1.0 - z*z)*dpy*dpy) ;
    w[n-i] = w[i-1] ;
  }

  return ;

}

