#include <Eigen/Dense>
#include <libint2.h>

#ifndef BASIS_H
#define BASIS_H

 std::vector<libint2::Shell> build_basis ( std::string & bas_name, Eigen::VectorXi AtN, Eigen::MatrixXd coord) ;

 std::vector<libint2::Shell> load_sto3g ( Eigen::VectorXi AtN, Eigen::MatrixXd c) ;

#endif
