#include <vector>

#ifndef INTEGR_H
#define INTEGR_H

void bonnet_r ( int n, float x, float& p2) ;

void bonnet_r( int n, float x, float& p2, float& dp) ;

void gauleg ( float lv, float uv, int n, std::vector<float>& x, std::vector<float>& w) ;

#endif
