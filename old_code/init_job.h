#include "common.h"
#include <Eigen/Core>
#include "nbodyint.h"
#include "tei.h"
#include <vector>

#ifndef INIT_JOB_H
#define INIT_JOB_H

template < class matrix>
void initialize( int ws, int hm, common& com, matrix& h, nbodyint<matrix>*& W, std::vector<tei>*& intarr, int& nbas) ;

#endif
