#include "basis.h"
#include "constants.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "tei.h"

#ifndef UTIL_H
#define UTIL_H

void bonnet_r ( int n, double x, double& p2) ;

void bonnet_r( int n, double x, double& p2, double& dp) ;

void charge_magnetic_decomposition( const Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> mx, Eigen::Ref<Eigen::MatrixXcd> my, Eigen::Ref<Eigen::MatrixXcd> mz) ;

double fact ( int n=1) ;

double factfact ( int n=1) ;

double fboys( int k, double t) ;

void fc_hamiltonian ( Eigen::Ref<Eigen::MatrixXcd> h, basis_set b, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::MatrixXd> n) ;

void transform( int o, Eigen::Ref<Eigen::MatrixXd> X, Eigen::Ref<Eigen::MatrixXd> M) ;

void transform( int o, Eigen::Ref<Eigen::MatrixXcd> X, Eigen::Ref<Eigen::MatrixXcd> M) ;

//void oao( Eigen::Ref<Eigen::MatrixXd> tmp, Eigen::Ref<Eigen::MatrixXd> O, Eigen::Ref<Eigen::MatrixXd> X) ;

void oao( int n, std::vector<tei>& iarr, std::vector<tei>& ioarr, Eigen::MatrixXd s) ;

void eulrgrd ( int n_ps, int n_t, int n_ph, std::vector<double>& w_ps, std::vector<double>& w_t, 
     std::vector<double>& w_ph, std::vector<double>& x_ps, std::vector<double>& x_t, 
     std::vector<double>& x_ph, int S) ;

//void K_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

//void F_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

//void T_op( hfwfn& a, Eigen::Ref<Eigen::MatrixXcd> m, int b ) ;

double pfaffian_A( Eigen::Ref<Eigen::MatrixXd> m) ;

cd pfaffian_A( Eigen::Ref<Eigen::MatrixXcd> m) ;

double pfaffian_H( Eigen::Ref<Eigen::MatrixXd> m) ;

cd pfaffian_H( Eigen::Ref<Eigen::MatrixXcd> m) ;

void rand_unitary( Eigen::Ref<Eigen::MatrixXcd> u) ;

void testing_magnetic_structure( const Eigen::Ref<Eigen::MatrixXcd> mx, const Eigen::Ref<Eigen::MatrixXcd> my, const Eigen::Ref<Eigen::MatrixXcd> mz) ;

void clean_mat( Eigen::Ref<Eigen::MatrixXcd> m) ;

uint32_t xorshift32( uint32_t *state) ;

double rand01( uint32_t *state) ;

void wavefunction_characterization( const int& nb,  const int& ne, const std::string& f = "qtz.wfn.bin") ;

#endif
