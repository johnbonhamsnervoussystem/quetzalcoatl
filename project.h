#include "common.h"
#include "diis.h"
#include "integr.h"
#include "nbodyint.h"
#include "tei.h"
#include <Eigen/Dense>
#include <vector>
#ifndef PROJECT_H
#define PROJECT_H

void prj_drv( common& com) ;

void proj_HFB( common& com) ;

//void rrHFB_projection( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, double& m, trapezoid*& ngrid, int& maxit, Eigen::Ref<Eigen::MatrixXcd> xs, Eigen::Ref<Eigen::MatrixXcd> xsi) ;

template < class matrix>
void projected_wavefunction_energy( const int& ne, const int& nb, nbodyint<matrix>* W, trapezoid*& ng, const Eigen::Ref<Eigen::MatrixXcd> h, Eigen::MatrixXcd rh, Eigen::MatrixXcd k, const Eigen::Ref<Eigen::MatrixXcd> k_i, Eigen::Ref<Eigen::MatrixXcd> Rp, Eigen::Ref<Eigen::MatrixXcd> rp, Eigen::Ref<Eigen::MatrixXcd> kp, Eigen::Ref<Eigen::MatrixXcd> kbp, Eigen::Ref<Eigen::MatrixXcd> I, Eigen::Ref<Eigen::MatrixXcd> s1, Eigen::Ref<Eigen::MatrixXcd> s2, Eigen::Ref<Eigen::MatrixXcd> s3, cd& X, cd& E) ;

template < class matrix>
void ring_shiekh_cg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, cd& norm, diis_control& diis_opt, int& maxit) ;

template < class matrix>
void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_opt, int& maxit) ;

template < class matrix>
void ring_shiekh_rr_K( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> V, Eigen::Ref<Eigen::MatrixXcd> U, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) ;

template < class matrix>
void Xring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_opt, int& maxit) ;

//void cgHFB_projection_dbg( int& ns, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& I, Eigen::Ref<Eigen::MatrixXcd> w, Eigen::Ref<Eigen::MatrixXcd> r, Eigen::Ref<Eigen::MatrixXcd> k, trapezoid*& ng, int& mit) ;

//cd pf_overlap( int& n, cd& c, Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, Eigen::Ref<Eigen::MatrixXcd> R) ;

void chemical_potential ( const int& nele, const int& nbas, double& lambda, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_diag, Eigen::Ref<Eigen::MatrixXcd> R, Eigen::Ref<Eigen::MatrixXcd> mu, Eigen::Ref<Eigen::MatrixXcd> eigvec, const Eigen::Ref<Eigen::MatrixXcd> Heff, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> H) ;

void jorge_guess( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, double n, double& norm) ;

#endif
