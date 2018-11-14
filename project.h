#include "common.h"
#include "hfwfn.h"
#include "integr.h"
#include "tei.h"
#include <Eigen/Dense>
#include <vector>
#ifndef PROJECT_H
#define PROJECT_H

void prj_drv( common& com, std::vector<tei>& intarr, int o) ;

void proj_HFB_cplx( common& com, std::vector<tei>& I) ;

void cgHFB_projection( int& ns, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& I, Eigen::Ref<Eigen::MatrixXcd> w, trapezoid*& ng, int& mit) ;

#endif

