#include<vector>
#include<Eigen/Dense>
#include "tei.h"

#ifndef SOLVER_H
#define SOLVER_H

void canort( const Eigen::Ref<Eigen::MatrixXf> s, Eigen::Ref<Eigen::MatrixXcf> xs , int& dim) ;

void canort( const Eigen::Ref<Eigen::MatrixXcf> s, Eigen::Ref<Eigen::MatrixXcf> xs , int& dim) ;

void scfdia ( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> const x, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> den) ;

void ahm_exp( Eigen::Ref<Eigen::MatrixXcf> x, Eigen::Ref<Eigen::MatrixXcf> u, int dim, int opt = 0) ;

#endif

