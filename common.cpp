#include<Eigen/Dense>
#include "common.h"

/* This stores information which is needed in most routines. */
/* Set the data */
  void common::nbas( int n) { nbasis = n ; return ; }
  void common::nbsu( int n) { nbsuse = n ; return ; }
  void common::ntei( int n) { n2ei   = n ; return ; }
  void common::natm( int n) { natoms = n ; return ; }
  void common::nele( int n) { nel    = n ; return ; }
  void common::nalp( int n) { nal    = n ; return ; }
  void common::nbet( int n) { nbe    = n ; return ; }
  void common::nrep( float f) { nn = f ; return ; }
 
/* Set matrix elements */
  void common::setS ( Eigen::MatrixXf s_in, int nb) {
    s_c.resize( nb, nb) ;
    s_c = s_in ;
    return ;
  }

  void common::setH ( Eigen::MatrixXf h_in, int nb) {
    h_c.resize( nb, nb) ;
    h_c = h_in ;
    return ;
  }


/* Retrieve the data */
  int common::nbas( void) {return nbasis ;}
  int common::nbsu( void) {return nbsuse ;}
  int common::ntei( void) {return n2ei   ;}
  int common::natm( void) {return natoms ;}
  int common::nele( void) {return nel    ;}
  int common::nalp( void) {return nal    ;}
  int common::nbet( void) {return nbe    ;}
  float common::nrep( void) {return nn   ;}

/* Retrieve a matrix element */
  float common::s( int i, int j) { return s_c( i, j) ;}
  float common::h( int i, int j) { return h_c( i, j) ;}

/* Retrieve a matrix */
  Eigen::MatrixXf common::getS( void) { return s_c ;}
  Eigen::MatrixXf common::getH( void) { return h_c ;}

