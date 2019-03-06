#include "constants.h"
#include "util.h"
#include <vector>

#ifndef INTEGR_H
#define INTEGR_H

class integration_grid {
  protected :
/*
  Weights and steps in grid
*/
    int steps ;
    std::vector<double> points ;
    std::vector<double> weights ;
    std::vector<double>::iterator ip ;
    std::vector<double>::iterator iw ;

  public :
/*
  Constructor
*/
  integration_grid( int n) {
    /*
      Allocate space for n points plus n = 0
    */
    steps = n ;
    points.reserve(n) ;
    weights.reserve(n) ;
    }
/*
  Access while iterating through a loop
*/
    int ns(void){ return steps ;}
    double q( void){ return *ip++ ;}
    double w( void){ return *iw++ ;}
/*
  Set the pointers to the beginning
*/
    void set_s( void){
      ip = points.begin() ;
      iw = weights.begin() ;
      return ;
      }
} ;

class trapezoid : public integration_grid {
  
  public :
    trapezoid ( double lv, double uv, int n) : integration_grid( n){
/* 
  Set up a trapezoid integration grid
*/
      double pt, wt ;
      double seg = ( uv - lv)/static_cast<double>(n) ;
  
      for ( auto i = 0; i <= n; i++){
        pt = seg*static_cast<double>(i) + lv ;
        weights.push_back(seg) ;
        points.push_back(pt) ;
        }
//      weights.push_back(d0) ;
//      points.push_back(d0) ;
      set_s() ;
      }
} ;


class gauleg : public integration_grid {

  public :
    gauleg ( double lv, double uv, int n) : integration_grid( n) {
/*
   Teukolsky, Vetterling and Flannery, Numerical Recipes in F77
   2nd Ed. 2005
 
   lv - lower value of the integration grid
   uv - upper value of the integration grid
   n - number of grid points
   x - abscissas of the quadrature
   w - weights of the quadrature
*/

      int m ;
      double z ;
      double z1 ;
      double xm ;
      double xl ;
      double py ;
      double dpy ;
      double eps = 1.0e-7 ;
      points.assign( n+1, d0) ;
      weights.assign( n+1, d0) ;

      m = (n + 1)/2 ;
      xm = (uv + lv)/2.0 ;
      xl = (uv - lv)/2.0 ;

      for ( int i=1; i <= m; i++ ) {

        z = std::cos ( pi*(static_cast<double>(i)-0.25)/(static_cast<double>(n)+0.5)) ;

        z1 = z + 1.0 ;
        dpy = 1.0 ;

        while ( std::abs(z-z1) > eps ) {

          bonnet_r ( n, z, py, dpy) ;
          z1 = z ;

          z = z1 - py/dpy ;
        }

        points[i-1] = xm - xl*z ;
        points[n-i] = xm + xl*z ;

        weights[i-1] = 2.0*xl/((1.0 - z*z)*dpy*dpy) ;
        weights[n-i] = weights[i-1] ;
      }
//      weights.push_back(d0) ;
//      points.push_back(d0) ;
      set_s() ;
    }
  } ; 

#endif
