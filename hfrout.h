#include "common.h"
#include<Eigen/Dense>
#include<vector>
#include<string>
#include "hfwfn.h"
#include "tei.h"

#ifndef HFROUT_H
#define HFROUT_H

void scf_drv( common& com, Eigen::Ref<Eigen::MatrixXd> T, Eigen::Ref<Eigen::MatrixXd> V, Eigen::Ref<Eigen::MatrixXd> S, std::vector<tei>& intarr, int& o);

double rrhfdia( Eigen::Ref<Eigen::MatrixXd> h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, int nbasis, int nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig ) ;

double crhfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig) ;

double ruhfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXd> c_a, Eigen::Ref<Eigen::MatrixXd> c_b, Eigen::Ref<Eigen::VectorXd> eig) ;

double cuhfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXcd> c_a, Eigen::Ref<Eigen::MatrixXcd> c_b, Eigen::Ref<Eigen::VectorXd> eig) ;

double rghfdia( Eigen::Ref<Eigen::MatrixXd> const h, Eigen::Ref<Eigen::MatrixXd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXd> c, Eigen::Ref<Eigen::VectorXd> eig ) ;

double cghfdia( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig) ;

#endif 
