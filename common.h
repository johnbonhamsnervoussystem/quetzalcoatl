#include <Eigen/Dense>
#include <vector>
#include <string>

#ifndef COMMON_H
#define COMMON_H

class common {
/* Job dimensions
 *
   store variables which are useful in most routines and will often 
   be used for dimensioning of matrices.

   nbasis - number of ao basis functions
   nbsuse - number of basis functions after orthogonalization of the 
            aos.
   ntei - the nuber of two electron integrals not related by symmetry
   natoms - the number of atoms.
   nel - the number of electrons
   nalp - the number of spin up electrons
   nbet - the number of spin down electrons

 * Stored Matrices
 *
   coord - coordinates matrix
   s - overlap matrix
   h - core hamiltonian

 */

private :
  int nbasis ;
  int n2ei   ;
  int natoms ;
  int nel    ;
  int nal    ;
  int nbe    ;
  double nn   ;
  std::string basis_name ;

  Eigen::MatrixXd s_c ;
  Eigen::MatrixXd h_c ;
  Eigen::VectorXi a_c ;
  Eigen::MatrixXd coord ;

public :
/* Set the data */
  void nbas( int n) ;
  void ntei( int n) ;
  void natm( int n) ;
  void nele( int n) ;
  void nalp( int n) ;
  void nbet( int n) ;
  void nrep( double f) ;
  void bnam( std::string n) ;

/* Coordinates */
  void setcoord( std::vector<std::vector<double>> c) ;

/* Set matrix elements */
  void setS( Eigen::MatrixXd s) ;
  void setH( Eigen::MatrixXd h) ;
  void setA( std::vector<int> a) ;
  void setC( std::vector<std::vector<double>> c) ;

/* Retrieve the data */
  int nbas( void) ;
  int ntei( void) ;
  int natm( void) ;
  int nele( void) ;
  int nalp( void) ;
  int nbet( void) ;
  double nrep( void) ;
  std::string bnam( void) ;

/* Retrieve a matrix */
  Eigen::MatrixXd getS( void) ;
  Eigen::MatrixXd getH( void) ;
  Eigen::VectorXi getA( void) ;
  Eigen::MatrixXd getC( void) ;

} ;
#endif
