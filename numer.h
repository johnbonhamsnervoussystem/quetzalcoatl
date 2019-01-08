#include <Eigen/Dense>
#include <vector>

#ifndef NUMER_H
#define NUMER_H

void numerical_derivative( double& h, Eigen::Ref<Eigen::MatrixXd> S, double (*f)( int&, cd&, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>), int& nbas, cd& c, Eigen::Ref<Eigen::MatrixXd> T, Eigen::Ref<Eigen::MatrixXd> U) ;

double five_point( const double& h, std::vector<double>& f) ;

#endif
