#include "constants.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "hfwfn.h"
#include "tei.h"

#ifndef UTIL_H
#define UTIL_H

void bonnet_r ( int n, double x, double& p2) ;

void bonnet_r( int n, double x, double& p2, double& dp) ;

double fact ( int n=1) ;

double factfact ( int n=1) ;

double fboys( int k, double t) ;

//template<class matrix>
//void oao( int& n, int& t, Eigen::MatrixBase<matrix>& a, Eigen::Ref<Eigen::MatrixXd> s, Eigen::Ref<Eigen::MatrixXd> x) ;

//void oao( Eigen::Ref<Eigen::MatrixXd> tmp, Eigen::Ref<Eigen::MatrixXd> O, Eigen::Ref<Eigen::MatrixXd> X) ;

void oao( int n, std::vector<tei>& iarr, std::vector<tei>& ioarr, Eigen::MatrixXd s) ;

void eulrgrd ( int n_ps, int n_t, int n_ph, std::vector<double>& w_ps, std::vector<double>& w_t, 
     std::vector<double>& w_ph, std::vector<double>& x_ps, std::vector<double>& x_t, 
     std::vector<double>& x_ph, int S) ;

void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

double pfaffian( Eigen::Ref<Eigen::MatrixXd> m) ;

cd pfaffian( Eigen::Ref<Eigen::MatrixXcd> m) ;

void rand_unitary( Eigen::Ref<Eigen::MatrixXcd> u) ;

#endif
