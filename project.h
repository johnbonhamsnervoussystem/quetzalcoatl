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

template < class matrix>
void general_derivative_testing( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::VectorXd> eigval, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) ;

template < class matrix>
void projected_wavefunction_energy( const int& ne, const int& nb, nbodyint<matrix>* W, trapezoid*& ng, const Eigen::Ref<Eigen::MatrixXcd> h, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kap, const Eigen::Ref<Eigen::MatrixXcd> kap_i, Eigen::Ref<Eigen::MatrixXcd> Rp, matrix& rp, matrix& kp, Eigen::Ref<Eigen::MatrixXcd> kbp, Eigen::Ref<Eigen::MatrixXcd> I, Eigen::Ref<Eigen::MatrixXcd> s1, Eigen::Ref<Eigen::MatrixXcd> s2, Eigen::Ref<Eigen::MatrixXcd> s3, cd& x_g, cd& E_g) ;

template <class matrix>
void build_projected_hamiltonian_energy( const int& nele, const int& nbas, nbodyint<matrix>* W, trapezoid*& ngrid, const Eigen::Ref<Eigen::MatrixXcd> h, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, const Eigen::Ref<Eigen::MatrixXcd> kappa_i, Eigen::Ref<Eigen::MatrixXcd> R_phi, matrix& r_phi, matrix& k_phi, Eigen::Ref<Eigen::MatrixXcd> kbar_phi, Eigen::Ref<Eigen::MatrixXcd> C_phi, Eigen::Ref<Eigen::MatrixXcd> G_phi, Eigen::Ref<Eigen::MatrixXcd> D_phi, Eigen::Ref<Eigen::MatrixXcd> I, Eigen::Ref<Eigen::MatrixXcd> epsilonN, Eigen::Ref<Eigen::MatrixXcd> gammaN, Eigen::Ref<Eigen::MatrixXcd> lambdaN, Eigen::Ref<Eigen::MatrixXcd> DeltaN, Eigen::Ref<Eigen::MatrixXcd> Y_tmp, Eigen::Ref<Eigen::MatrixXcd> Y_phi, Eigen::Ref<Eigen::MatrixXcd> Heff, Eigen::Ref<Eigen::MatrixXcd> t, Eigen::Ref<Eigen::MatrixXcd> tmp, cd& intx_g, cd& intE_g) ;

template < class matrix>
void ring_shiekh_cg( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, cd& norm, diis_control& diis_opt, int& maxit) ;

template < class matrix>
void ring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::VectorXd> eval, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_opt, int& maxit) ;

template < class matrix>
void ring_shiekh_K( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> V, Eigen::Ref<Eigen::MatrixXcd> U, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_cntrl, int& maxit) ;

template < class matrix>
void Xring_shiekh_rr( int& nbas, Eigen::Ref<Eigen::MatrixXcd> h, nbodyint<matrix>* W, Eigen::Ref<Eigen::MatrixXcd> c, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> kappa, double& lambda, cd lshift, trapezoid*& ngrid, double& nele, diis_control& diis_opt, int& maxit) ;

void chemical_potential ( const int& nele, const int& nbas, double& lambda, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>& H_diag, Eigen::Ref<Eigen::MatrixXcd> R, Eigen::Ref<Eigen::MatrixXcd> mu, Eigen::Ref<Eigen::MatrixXcd> eigvec, const Eigen::Ref<Eigen::MatrixXcd> Heff, Eigen::Ref<Eigen::MatrixXcd> rho, Eigen::Ref<Eigen::MatrixXcd> H) ;

void jorge_guess( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> k, double n, double& norm) ;

#endif
