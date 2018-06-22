#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <vector>
#include "common.h"
#include "hfwfn.h"
#include "tei.h"

#ifndef SOLVER_H
#define SOLVER_H

void canort( const Eigen::Ref<Eigen::MatrixXf> s, Eigen::Ref<Eigen::MatrixXcf> xs , int& dim) ;

void canort( const Eigen::Ref<Eigen::MatrixXcf> s, Eigen::Ref<Eigen::MatrixXcf> xs , int& dim) ;

void scfdia ( Eigen::Ref<Eigen::MatrixXcf> const h, Eigen::Ref<Eigen::MatrixXcf> const x, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> den) ;

void ahm_exp( Eigen::Ref<Eigen::MatrixXcf> x, Eigen::Ref<Eigen::MatrixXcf> u, int dim, int opt = 0) ;

void trci( common& c, std::vector<hfwfn>& d, Eigen::Ref<Eigen::MatrixXcf> H, std::vector<tei>& intarr) ;

void trci( common& com, std::vector<hfwfn>& det, Eigen::Ref<Eigen::MatrixXcf> H, std::vector<tei>& intarr, std::string& trd, std::string& fop) ;

#endif

