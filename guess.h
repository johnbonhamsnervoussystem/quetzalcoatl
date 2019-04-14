#include "common.h"
#include <Eigen/Dense>

#ifndef GUESS_H
#define GUESS_H

void generate_pk( int nb, int ne, Eigen::Ref<Eigen::MatrixXcd> r, Eigen::Ref<Eigen::MatrixXcd> k) ;

void guess_drv( common& c) ;

template< class matrix>
void homo_lumo_mix_r( matrix& w, const int i, double angle) ;

void homo_lumo_mix_c( Eigen::Ref<Eigen::MatrixXcd> w, const int i, double angle) ;

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k) ;

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> k) ;

void Uthermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k) ;

void rand01_guess( Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) ;

void thermal_occm( int& nb, int& ne, Eigen::Ref<Eigen::VectorXd> eig, double& T, Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) ;

void thermal_occv( int& nbas, int& nele, Eigen::Ref<Eigen::VectorXd> eig, double& T, Eigen::Ref<Eigen::VectorXd> A, Eigen::Ref<Eigen::VectorXd> B) ;

#endif
