#include "common.h"
#include "integr.h"
#include "tei.h"
#include <Eigen/Dense>
#include <vector>
#ifndef PROJECT_H
#define PROJECT_H

void prj_drv( common& com, std::vector<tei>& intarr, int o) ;

void proj_HFB_real( common& com, std::vector<tei>& I) ;

void proj_HFB_cplx( common& com, std::vector<tei>& I) ;

void rrHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& m, trapezoid*& ngrid, int& maxit, Eigen::Ref<Eigen::MatrixXcd> Xsc, Eigen::Ref<Eigen::MatrixXcd> Xsci) ;

void cgHFB_projection_dbg( int& ns, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& I, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> r, Eigen::Ref<Eigen::MatrixXcd> k, trapezoid*& ng, int& mit) ;

cd pf_overlap( int& n, cd& c, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, Eigen::Ref<Eigen::MatrixXcd> R) ;

#endif

