#include "common.h"
#include <Eigen/Dense>
#include <vector>
#include <string>

/* Intialize defaults */
  common::common( void) {
    max_scf_iter = 20 ;
    scf_convergence_threshold = 1.0e-9 ;
    return ;
  } 
/* This stores information which is needed in most routines.
 Set things --
  General data for calculations */
  void common::nbas( int n) { nbasis = n ; return ;}
  void common::ntei( int n) { n2ei   = n ; return ;}
  void common::natm( int n) { natoms = n ; return ;}
  void common::nele( int n) { nel    = n ; return ;}
  void common::nalp( int n) { nal    = n ; return ;}
  void common::nbet( int n) { nbe    = n ; return ;}
  void common::nrep( double f) { nn = f ; return ;}
  void common::bnam( std::string n) { basis_name = n ; return ; }

/* Routine/algorithm control options */
  void common::mxscfit( int n) { max_scf_iter = n ; return ;}
  void common::scfthresh( double d) { scf_convergence_threshold = d ; return ;}
  

/* Matrix elements */
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
  void common::setA ( std::vector<double> a) {
    /* Save the atomic coordinates into an eigen array */
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

/* Get things */
  int common::nbas( void) {return nbasis ;}
  int common::ntei( void) {return n2ei   ;}
  int common::natm( void) {return natoms ;}
  int common::nele( void) {return nel    ;}
  int common::nalp( void) {return nal    ;}
  int common::nbet( void) {return nbe    ;}
  double common::nrep( void) {return nn   ;}
  std::string common::bnam( void) {return basis_name ;}

  int common::mxscfit( void) {return max_scf_iter ;}
  double common::scfthresh( void) {return scf_convergence_threshold ;}

/* Retrieve a matrix */
  Eigen::MatrixXd common::getS( void) { return s_c   ;}
  Eigen::MatrixXd common::getH( void) { return h_c   ;}
  Eigen::VectorXd common::getA( void) { return a_c   ;}
  Eigen::MatrixXd common::getC( void) { return coord ;}

