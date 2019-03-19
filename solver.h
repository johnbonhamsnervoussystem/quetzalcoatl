#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <vector>
#include "common.h"
#include "tei.h"

#ifndef SOLVER_H
#define SOLVER_H

void symort( Eigen::Ref<Eigen::MatrixXd> const s, Eigen::Ref<Eigen::MatrixXcd> T) ;

void canort( const Eigen::Ref<Eigen::MatrixXd> s, Eigen::Ref<Eigen::MatrixXcd> xs , int dim) ;

void canort( const Eigen::Ref<Eigen::MatrixXcd> s, Eigen::Ref<Eigen::MatrixXcd> xs , int dim) ;

void scfdia ( Eigen::Ref<Eigen::MatrixXcd> const h, Eigen::Ref<Eigen::MatrixXcd> const x, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> den) ;

void ahm_exp( Eigen::Ref<Eigen::MatrixXcd> x, Eigen::Ref<Eigen::MatrixXcd> u, int dim, int opt = 0) ;

//void trci( common& c, std::vector<hfwfn>& d, Eigen::Ref<Eigen::MatrixXcd> H, std::vector<tei>& intarr) ;

//void trci( common& com, std::vector<hfwfn>& det, Eigen::Ref<Eigen::MatrixXcd> H, std::vector<tei>& intarr, std::string& trd, std::string& fop) ;

#endif

