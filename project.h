#include "common.h"
#include "integr.h"
#include "nbodyint.h"
#include "tei.h"
#include <Eigen/Dense>
#include <vector>
#ifndef PROJECT_H
#define PROJECT_H

void prj_drv( common& com) ;

void proj_HFB_real( common& com) ;

void proj_HFB_cplx( common& com, std::vector<tei>& I) ;

//void rrHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, double& m, trapezoid*& ngrid, int& maxit, Eigen::Ref<Eigen::MatrixXcd> xs, Eigen::Ref<Eigen::MatrixXcd> xsi) ;

template < class matrix>
void ring_shiekh_cg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, cd& norm, int& maxit) ;

template < class matrix>
void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::VectorXd> eig, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, int& maxit) ;

template < class matrix>
void tr_ring_shiekh( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, int& maxit) ;

//void cgHFB_projection_dbg( int& ns, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& I, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> r, Eigen::Ref<Eigen::MatrixXcd> k, trapezoid*& ng, int& mit) ;

//cd pf_overlap( int& n, cd& c, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, Eigen::Ref<Eigen::MatrixXcd> R) ;

void jorge_guess( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, double n, double& norm) ;

#endif

