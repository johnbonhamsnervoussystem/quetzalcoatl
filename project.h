#include "common.h"
#include "hfwfn.h"
#include "tei.h"
#include <Eigen/Dense>
#include <vector>
#ifndef PROJECT_H
#define PROJECT_H

void prj_drv( common& com, int o) ;

void proj_HFB_cplx( common& com) ;

void e_phf( common& com, hfwfn& refd, Eigen::Ref<Eigen::MatrixXcd> H, std::vector<tei>& intarr) ;

#endif

