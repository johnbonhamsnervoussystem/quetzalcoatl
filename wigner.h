#include "constants.h"
#include "common.h"
#include <Eigen/Dense>

#ifndef WIGNER_H
#define WIGNER_H

double d_cof ( int j, int m, int k, int n) ;

double d_cof ( double j, double m, double k, int n) ;

double small_wd ( int j, int m, int k, double beta) ;

double small_wd ( double j, double m, double k, double beta) ;

std::complex<double> wigner_D( int j, int m, int k, double alpha, double beta, double gamma) ;
 
std::complex<double> wigner_D( double j, double m, double k, double alpha, double beta, double gamma) ;

//void R_s ( common& c, hfwfn& a, hfwfn& b, double al, double be, double ga) ;
//
//void R_s ( int n, hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, double al, double be, double ga) ;

void R_s ( int n, Eigen::Ref<Eigen::MatrixXcd> moa, Eigen::Ref<Eigen::MatrixXcd> mob, double al, double be, double ga) ;

#endif
