#include "common.h"
#include <Eigen/Dense>
#include "hfwfn.h"
#include <string>
#include "tei.h"
#include <vector>
#ifndef QTZIO_H
#define QTZIO_H

void read_input( common& c, const std::string& in_s) ;

void getmel( std::string file1, std::string file2, std::vector<tei>& I, common& c) ;

void rdsdet ( int n, std::vector<std::string>& m, std::vector<hfwfn>& d) ;

#endif
