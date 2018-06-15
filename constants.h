#include <complex>
#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Type definitions and values which are generally useful. */
/* Here is the wishlist but for sanity, tackle this one constant
 * at a time.
 * 
 * typedef cf complex<float>
 * d0 = 0.0 
 * d1 = 1.0
 * d2 = 2.0
 * zi = complex( 0.0, 1.0)
 * */

  typedef std::complex<float> cf ;
  extern const float pi ;
  extern const float d0 ;
  extern const cf z0 ;
  extern const cf z1 ;
  extern const cf zi ;

#endif
