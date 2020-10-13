#include <Eigen/Core>
#include <fstream>
#include <string>

#ifndef BINIO_H
#define BINIO_H

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ostream& F_OUT) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::istream& F_IN) ;

#endif

