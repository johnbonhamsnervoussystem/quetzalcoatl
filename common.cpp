#include "common.h"
#include <Eigen/Dense>
#include <vector>
#include <string>

/* This stores information which is needed in most routines. */
/* Set the data */
  void common::nbas( int n) { nbasis = n ; return ;}
  void common::ntei( int n) { n2ei   = n ; return ;}
  void common::natm( int n) { natoms = n ; return ;}
  void common::nele( int n) { nel    = n ; return ;}
  void common::nalp( int n) { nal    = n ; return ;}
  void common::nbet( int n) { nbe    = n ; return ;}
  void common::nrep( double f) { nn = f ; return ;}
  void common::bnam( std::string n) { basis_name = n ; return ; }

/* Set matrix elements */
  void common::setS ( Eigen::MatrixXd s_in) {
    s_c.resize( nbasis, nbasis) ;
    s_c = s_in ;
    return ;
  }

  void common::setH ( Eigen::MatrixXd h_in) {
    h_c.resize( nbasis, nbasis) ;
    h_c = h_in ;
    return ;
  }

/* Atomic Number */
  void common::setA ( std::vector<int> a) {
    /* Save the coordinates into an Eigen Array */
    a_c.resize(natoms) ;
    for ( int i=0; i < a.size(); i++){
      a_c[i] = a[i] ;
    }
    return ;
  }

/* Coordinates */
  void common::setC ( std::vector<std::vector<double>> c) {
    /* Save the coordinates into an Eigen Array */
    coord.resize(natoms,3) ;
    for ( int i=0; i < c.size(); i++){
      coord.row(i) << c[i][0], c[i][1], c[i][2] ;
    }
    return ;
  }

/* Retrieve the data */
  int common::nbas( void) {return nbasis ;}
  int common::ntei( void) {return n2ei   ;}
  int common::natm( void) {return natoms ;}
  int common::nele( void) {return nel    ;}
  int common::nalp( void) {return nal    ;}
  int common::nbet( void) {return nbe    ;}
  double common::nrep( void) {return nn   ;}
  std::string common::bnam( void) {return basis_name ;}

/* Retrieve a matrix */
  Eigen::MatrixXd common::getS( void) { return s_c   ;}
  Eigen::MatrixXd common::getH( void) { return h_c   ;}
  Eigen::VectorXi common::getA( void) { return a_c   ;}
  Eigen::MatrixXd common::getC( void) { return coord ;}

