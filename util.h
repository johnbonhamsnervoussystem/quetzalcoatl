#include<Eigen/Dense>
#include<vector>
#include "hfwfn.h"
#include "tei.h"

#ifndef UTIL_H
#define UTIL_H

float fact ( int n=1) ;

void oao( int n, hfwfn& a, Eigen::MatrixXf s) ;

void oao( int n, Eigen::Ref<Eigen::MatrixXf> ouv, Eigen::MatrixXf s) ;

void oao( int n, std::vector<tei>& iarr, std::vector<tei>& ioarr, Eigen::MatrixXf s) ;

void eulrgrd ( int n_ps, int n_t, int n_ph, std::vector<float>& w_ps, std::vector<float>& w_t, 
     std::vector<float>& w_ph, std::vector<float>& x_ps, std::vector<float>& x_t, 
     std::vector<float>& x_ph, int S) ;

void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcf> m, int b ) ;

void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcf> m, int b ) ;

void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcf> m, int b ) ;

#endif
