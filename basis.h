#include <Eigen/Dense>
#include <vector>
#include <string>

#ifndef BASIS_H
#define BASIS_H

  struct gau_prm {
    /* Gaussian primitive data
      x - Exponent
      c - coefficient
    */
    double x ;
    double c ;
    } ;

  struct sto {
    /* Slater Type Orbital
       l - angular momentum of x,y,z
       nprm - number of primitives
       g - vector of primitives
    */
    Eigen::Vector3i l ;
    int nprm ;
    double norm ;
    std::vector<gau_prm> g ;
    } ;

  struct basis_fc {
    /* Basis function
       nshl - number of shells( number of Stos)
       c - coordinates
       s - vector of sto*/
    int nshl ;
    Eigen::Vector3d c ;
    std::vector<sto> s ;
    } ;

  struct basis_set {
    /* Container for all the basis functions */
    int nbas ;
    std::vector<basis_fc> b ;
    } ;

  basis_set load_sto3g ( Eigen::VectorXi A, Eigen::MatrixXd c) ;

  basis_set build_basis ( std::string n, Eigen::Ref<Eigen::VectorXi> A, Eigen::Ref<Eigen::MatrixXd> c) ;

  void norm_basis ( int n, basis_set& b) ;

  double norm_gprm ( double a, Eigen::Vector3i L) ;

  void print_basis( int n, basis_set& b) ;

#endif
