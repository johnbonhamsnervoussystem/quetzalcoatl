#include "common.h"
#include <Eigen/Core>
#include <fstream>
#include "hfwfn.h"
#include <string>
#include "tei.h"
#include <vector>
#ifndef QTZIO_H
#define QTZIO_H

void read_input( common& c, const std::string& in_s) ;

void getmel( std::string file1, std::string file2, std::vector<tei>& I, common& c) ;

void strip_lower( std::string& s) ;

//void rdsdet ( int n, std::vector<std::string>& m, std::vector<hfwfn>& d) ;

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ofstream& F_OUT) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::ifstream& F_IN) ;

#endif
