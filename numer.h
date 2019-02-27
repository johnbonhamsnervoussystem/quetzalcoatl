#include <Eigen/Dense>
#include <vector>

#ifndef NUMER_H
#define NUMER_H

void numerical_derivative( double& h, Eigen::Ref<Eigen::MatrixXd> S, double (*f)( int&, cd&, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>), int& n, cd& c, Eigen::Ref<Eigen::MatrixXd> T, Eigen::Ref<Eigen::MatrixXd> U) ;

void numerical_derivative( double& h, Eigen::Ref<Eigen::MatrixXcd> S, cd (*f)( int&, cd&, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>, Eigen::Ref<Eigen::MatrixXcd>), int& n, cd& c, Eigen::Ref<Eigen::MatrixXcd> T, Eigen::Ref<Eigen::MatrixXcd> U) ;

double five_point( const double& h, std::vector<double>& f) ;

cd five_point( const double& h, std::vector<cd>& f) ;

#endif
