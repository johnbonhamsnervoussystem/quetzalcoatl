#include<Eigen/Dense>
#include<vector>
#include<string>
#include "hfwfn.h"
#include "tei.h"

#ifndef HFROUT_H
#define HFROUT_H

float rrhfdia( Eigen::Ref<Eigen::MatrixXf> h, Eigen::Ref<Eigen::MatrixXf> s, std::vector<tei>& intarr, int nbasis, int nele, Eigen::Ref<Eigen::MatrixXf> c, Eigen::Ref<Eigen::VectorXf> eig ) ;

float crhfdia( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcf> c, Eigen::Ref<Eigen::VectorXf> eig) ;

float ruhfdia( Eigen::Ref<Eigen::MatrixXf> const h, Eigen::Ref<Eigen::MatrixXf> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXf> c_a, Eigen::Ref<Eigen::MatrixXf> c_b, Eigen::Ref<Eigen::VectorXf> eig) ;

float cuhfdia( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> s, std::vector<tei>& intarr, const int& nbasis, const int& nalp, const int& nbet, Eigen::Ref<Eigen::MatrixXcf> c_a, Eigen::Ref<Eigen::MatrixXcf> c_b, Eigen::Ref<Eigen::VectorXf> eig) ;

float rghfdia( Eigen::Ref<Eigen::MatrixXf> const h, Eigen::Ref<Eigen::MatrixXf> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXf> c, Eigen::Ref<Eigen::VectorXf> eig ) ;

float cghfdia( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> s, std::vector<tei>& intarr, const int& nbasis, const int& nele, Eigen::Ref<Eigen::MatrixXcf> c, Eigen::Ref<Eigen::VectorXf> eig) ;

#endif 
