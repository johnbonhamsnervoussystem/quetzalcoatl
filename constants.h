#include <complex>
#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Type definitions and values which are generally useful. */
/* Here is the wishlist but for sanity, tackle this one constant
 * at a time.
 * 
 * d1 = 1.0
 * d2 = 2.0
 * zi = complex( 0.0, 1.0)
 * */

  typedef std::complex<double> cd ;
  extern const double pi ;
  extern const double d0 ;
  extern const double d1 ;
  extern const double d2 ;
  extern const cd z0 ;
  extern const cd z1 ;
  extern const cd zi ;

#endif
