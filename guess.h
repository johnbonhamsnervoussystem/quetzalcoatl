#include "common.h"
#include <Eigen/Dense>

#ifndef GUESS_H
#define GUESS_H

void guess_drv( common& c) ;

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k) ;

void thermal_guess( int& nele, int& nbas, Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> k) ;

void rand01_guess( Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) ;

void thermal_occ( int& nb, int& ne, Eigen::Ref<Eigen::VectorXd> eig, double& T, Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::MatrixXd> B) ;

#endif
