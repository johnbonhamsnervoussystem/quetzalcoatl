#include "constants.h"
#include <cmath>
#include <complex>

  /* Useful computing quantities */
  const double d0 = 0.0e0 ;
  const double d1 = 1.0e0 ;
  const double d2 = 2.0e0 ;
  const double d3 = 3.0e0 ;
  const double d4 = 4.0e0 ;
  const double d5 = 5.0e0 ;
  const double d6 = 6.0e0 ;
  const double d7 = 7.0e0 ;
  const double d8 = 8.0e0 ;
  const double d9 = 9.0e0 ;
  const double d10 = 10.0e0 ;
  const cd zpi ( pi, d0) ;
  const cd z0 = cd( 0.0e0, 0.0e0) ;
  const cd z1 = cd( 1.0e0, 0.0e0) ;
  const cd z2 = cd( 2.0e0, 0.0e0) ;
  const cd z3 = cd( 3.0e0, 0.0e0) ;
  const cd z4 = cd( 4.0e0, 0.0e0) ;
  const cd z5 = cd( 5.0e0, 0.0e0) ;
  const cd z6 = cd( 6.0e0, 0.0e0) ;
  const cd z7 = cd( 6.0e0, 0.0e0) ;
  const cd z8 = cd( 8.0e0, 0.0e0) ;
  const cd z9 = cd( 9.0e0, 0.0e0) ;
  const cd z10 = cd( 10.0e0, 0.0e0) ;
  const cd zi = cd( 0.0e0, 1.0e0) ;

  /* Useful physical quantities */
  extern const double bh2an = 0.5291772086 ;
  /* Boltzmann's Constatn in Hartree/Kelvin */
  extern const double kb = 3.1668114e-6 ;
  /* Pi */
  const double pi = d4*std::atan(d1) ;
