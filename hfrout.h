#include "common.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "nbodyint.h"
#include "tei.h"

#ifndef HFROUT_H
#define HFROUT_H

void scf_drv( common& c) ;

void real_HFB( common& c, int& o) ;

void cplx_HFB( common& c, int& o) ;

void real_SlaDet( common& c, int& o) ;

void cplx_SlaDet( common& c, int& o) ;

template < class matrix>
double rhfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nele, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, int diistype, bool ls, double lvs, const int& maxit, const double& thresh) ;

template < class matrix>
double uhfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nalp, const int& nbet, matrix& c_a, matrix& c_b, Eigen::Ref<Eigen::VectorXd> eig, const int& maxit, const double& thresh) ;

template < class matrix>
double ghfdia( const matrix& h, nbodyint<matrix>* W, const int& nbasis, const int& nele, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, int n, bool ls, double lvs, const int& maxit, const double& thresh) ;

double ghfdia_fc( const Eigen::Ref<Eigen::MatrixXcd> h, const Eigen::Ref<Eigen::MatrixXcd> f, nbodyint<Eigen::MatrixXcd>* W, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, int n, bool lshift, double lvlshft, const int& maxit, const double& thresh) ;

template < class matrix, class z>
double rhfbdia( const matrix& h, nbodyint<matrix>* X, const int& nbasis, const int& nele, matrix& p, matrix& k, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, z& lambda, z& lshift, const int& maxit_scf, const int& maxit_pn, const double& thresh) ;

template < class matrix, class z>
double ghfbdia( const matrix& h, nbodyint<matrix>* X, const int& nbasis, const int& nele, matrix& p, matrix& k, matrix& c, Eigen::Ref<Eigen::VectorXd> eig, z& lambda, z& lshift, const int& maxit_scf, const int& maxit_pn, const double& thresh) ;

#endif
