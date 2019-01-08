#include "constants.h"
#include <cmath>
#include <complex>

  /* Useful computing quantities */
  const double pi = 4.0e0*std::atan(1.0e0) ;
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
  const cd z0 = cd( 0.0e0, 0.0e0) ;
  const cd z1 = cd( 1.0e0, 0.0e0) ;
  const cd z2 = cd( 2.0e0, 0.0e0) ;
  const cd zi = cd( 0.0e0, 1.0e0) ;

  /* Useful physical quantities */
  extern const double bh2an = 0.5291772086 ;
  /* Boltzmann's Constatn in Hartree/Kelvin */
  extern const double kb = 3.1668114e-6 ;
