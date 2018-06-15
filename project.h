#include "common.h"
#include "hfwfn.h"
#include "tei.h"
#include <Eigen/Dense>
#include <vector>
#ifndef PROJECT_H
#define PROJECT_H

void e_phf( common& com, hfwfn& refd, Eigen::Ref<Eigen::MatrixXcf> H, std::vector<tei>& intarr) ;

#endif

