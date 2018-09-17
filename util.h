#include<Eigen/Dense>
#include<vector>
#include "hfwfn.h"
#include "tei.h"

#ifndef UTIL_H
#define UTIL_H

double fact ( int n=1) ;

double factfact ( int n=1) ;

double fboys( int k, double t) ;

void oao( int n, hfwfn& a, Eigen::MatrixXd s) ;

void oao( int n, Eigen::Ref<Eigen::MatrixXd> ouv, Eigen::MatrixXd s) ;

void oao( int n, std::vector<tei>& iarr, std::vector<tei>& ioarr, Eigen::MatrixXd s) ;

void eulrgrd ( int n_ps, int n_t, int n_ph, std::vector<double>& w_ps, std::vector<double>& w_t, 
     std::vector<double>& w_ph, std::vector<double>& x_ps, std::vector<double>& x_t, 
     std::vector<double>& x_ph, int S) ;

void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

double pfaffian( Eigen::Ref<Eigen::MatrixXd> m) ;

#endif
