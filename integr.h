#include <vector>

#ifndef INTEGR_H
#define INTEGR_H

void bonnet_r ( int n, double x, double& p2) ;

void bonnet_r( int n, double x, double& p2, double& dp) ;

void gauleg ( double lv, double uv, int n, std::vector<double>& x, std::vector<double>& w) ;

#endif
