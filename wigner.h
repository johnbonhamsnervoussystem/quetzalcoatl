#include "hfwfn.h"
#include "common.h"

#ifndef WIGNER_H
#define WIGNER_H

float d_cof ( int j, int m, int k, int n) ;

float d_cof ( float j, float m, float k, int n) ;

float small_wd ( int j, int m, int k, float beta) ;

float small_wd ( float j, float m, float k, float beta) ;

std::complex<float> wigner_D( int j, int m, int k, float alpha, float beta, float gamma) ;
 
std::complex<float> wigner_D( float j, float m, float k, float alpha, float beta, float gamma) ;

void R_s ( common& c, hfwfn& a, hfwfn& b, float alpha, float beta, float gamma) ;

#endif
