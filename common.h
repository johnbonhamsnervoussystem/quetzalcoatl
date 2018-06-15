#include<Eigen/Dense>

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
   s - overlap matrix
   h - core hamiltonian

 */

private :
  int nbasis ;
  int nbsuse ;
  int n2ei   ;
  int natoms ;
  int nel    ;
  int nal    ;
  int nbe    ;
  float nn   ;

  Eigen::MatrixXf s_c ;
  Eigen::MatrixXf h_c ;

public :
/* Set the data */
  void nbas( int n) ;
  void nbsu( int n) ;
  void ntei( int n) ;
  void natm( int n) ;
  void nele( int n) ;
  void nalp( int n) ;
  void nbet( int n) ;
  void nrep( float f) ;

/* Set matrix elements */
  void setS( Eigen::MatrixXf s_in, int nb) ;
  void setH( Eigen::MatrixXf h_in, int nb) ;

/* Retrieve the data */
  int nbas( void) ;
  int nbsu( void) ;
  int ntei( void) ;
  int natm( void) ;
  int nele( void) ;
  int nalp( void) ;
  int nbet( void) ;
  float nrep( void) ;

/* Access matrix elements */
  float s( int i, int j) ;
  float h( int i, int j) ;

/* Retrieve a matrix */
  Eigen::MatrixXf getS( void) ;
  Eigen::MatrixXf getH( void) ;

} ;
#endif
