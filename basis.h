#include <Eigen/Dense>

#ifndef BASIS_H
#define BASIS_H

 std::vector<libint2::Shell> build_basis ( std::string & bas_name, Eigen::VectorXi AtN, Eigen::MatrixXd coord) ;

#endif
