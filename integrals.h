#include <vector>
#include <libint2.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix; 

void compute_onebody(libint2::Engine e, libint2::BasisSet b, Eigen::MatrixXd& m);
