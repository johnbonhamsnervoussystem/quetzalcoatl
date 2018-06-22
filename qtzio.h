#include<string>
#include<vector>
#include<Eigen/Dense>
#include "common.h"
#include "hfwfn.h"
#include "tei.h"
#ifndef QTZIO_H
#define QTZIO_H

void getmel(std::string file1, std::string file2, std::vector<tei>& I, common& c) ;

void rdsdet ( int n, std::vector<std::string>& m, std::vector<hfwfn>& d) ;

#endif
