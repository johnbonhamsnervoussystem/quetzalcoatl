#include <Eigen/Core>
#include <fstream>
#include "binio.h"

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ostream& F_OUT) {
  typename matrix::Index rows=m.rows(), cols=m.cols();

  F_OUT.write( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;

  return ;

}

template void write_eigen_bin(const Eigen::MatrixXf&, std::ostream&) ;
template void write_eigen_bin(const Eigen::MatrixXd&, std::ostream&) ;
template void write_eigen_bin(const Eigen::MatrixXcf&, std::ostream&) ;
template void write_eigen_bin(const Eigen::MatrixXcd&, std::ostream&) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::istream& F_IN) {
  typename matrix::Index rows=m.rows(), cols=m.cols();

  F_IN.read( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;

  return ;

}

template void read_eigen_bin(const Eigen::MatrixXf&, std::istream&) ;
template void read_eigen_bin(const Eigen::MatrixXd&, std::istream&) ;
template void read_eigen_bin(const Eigen::MatrixXcf&, std::istream&) ;
template void read_eigen_bin(const Eigen::MatrixXcd&, std::istream&) ;

