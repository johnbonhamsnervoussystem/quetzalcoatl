#include "common.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "tei.h"

#ifndef HFROUT_H
#define HFROUT_H

void scf_drv( common& com, std::vector<tei>& intarr, int o) ;

void real_HFB( common& com, std::vector<tei>& intarr, int opt) ;

void cplx_HFB( common& com, std::vector<tei>& intarr, int opt) ;

void real_SlaDet( common& com, std::vector<tei>& intarr, int o) ;

void cplx_SlaDet( common& com, std::vector<tei>& intarr, int o) ;

double rrhfdia( Eigen::Ref<Eigen::MatrixXd> h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, int nbasis, int nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

double crhfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

double ruhfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXd> c_a, Eigen::Ref<Eigen::MatrixXd> c_b, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

double cuhfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXcd> c_a, Eigen::Ref<Eigen::MatrixXcd> c_b, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

double rghfbdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

double rghfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

double rrhfbdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> k, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig, double& lambda, int& maxit_scf, int& maxit_pn, double& thresh) ;

double cghfbdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, double& l, int& mis, int& min, double& t) ;

double cghfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, int& mi, double& t) ;

#endif 
