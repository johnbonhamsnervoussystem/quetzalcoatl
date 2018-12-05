#include <Eigen/Dense>
#include <vector>
#include <string>

#ifndef COMMON_H
#define COMMON_H

/*

  At some point I'd like not to intialize values I won't use.  For
  now things are small enough that it is likely not an issue.

*/

class common {
/* 
  Job Information
    hamiltonian - Indiacates the type of matrix elements we are computing
    method - What type of job are we doing.

  System information :
    nbasis - number of ao basis functions
    nbsuse - number of basis functions after orthogonalization of the 
             aos.
    ntei - the nuber of two electron integrals not related by symmetry
    natoms - the number of atoms.
    nel - the number of electrons
    nalp - the number of spin up electrons
    nbet - the number of spin down electrons
    chemical_potential - the initial guess for mu

  Algorithm Control :
    max_scf_iter - the maximum number of scf iterations
    scf_convergence_threshold - the threshold at which we consider the wavefunction
      converged.

  Stored Matrices
   s - overlap matrix
   h - core hamiltonian
   a - atomic number
   coord - coordinates matrix

  Projection Data
    ngrid - particle number projection grid
*/

private :

  int hamiltonian ;
  int method ;
  int nbasis ;
  int n2ei ;
  int natoms ;
  int nel ;
  int nal ;
  int nbe ;
  double chemical_potential ;
  double nn ;
  std::string basis_name ;

  int max_scf_iter ;
  double scf_convergence_threshold ;

  Eigen::MatrixXd s_c ;
  Eigen::MatrixXd xs_c ;
  Eigen::MatrixXd h_c ;
  Eigen::VectorXd a_c ;
  Eigen::MatrixXd coord ;

  int nbr_g ;

public :
/* Initializers to set default values */
  common( void) ;

/* Set the data */
  void hamil( int n) ;
  void methd( int n) ;
  void nbas( int n) ;
  void ntei( int n) ;
  void natm( int n) ;
  void nele( int n) ;
  void nalp( int n) ;
  void nbet( int n) ;
  void mu( int n) ;
  void nrep( double f) ;
  void bnam( std::string n) ;
  void ngrid( int n) ;

/* Algorithm control*/
  void mxscfit( int n) ;
  void scfthresh( double d) ;

/* Coordinates */
  void setcoord( std::vector<std::vector<double>> c) ;

/* Set matrix elements */
  void setS( Eigen::MatrixXd s) ;
  void setXS( Eigen::MatrixXd xs) ;
  void setH( Eigen::MatrixXd h) ;
  void setA( std::vector<double> a) ;
  void setC( std::vector<std::vector<double>> c) ;

/* Retrieve the data */
  int hamil( void) ;
  int methd( void) ;
  int nbas( void) ;
  int ntei( void) ;
  int natm( void) ;
  int nele( void) ;
  int nalp( void) ;
  int nbet( void) ;
  double mu( void) ;
  double nrep( void) ;
  std::string bnam( void) ;
/* Projection Data */
  int ngrid( void) ;

/* Algorithm control*/
  int mxscfit( void) ;
  double scfthresh( void) ;

/* Retrieve a matrix */
  Eigen::MatrixXd getS( void) ;
  Eigen::MatrixXd getXS( void) ;
  Eigen::MatrixXd getH( void) ;
  Eigen::VectorXd getA( void) ;
  Eigen::MatrixXd getC( void) ;


} ;
#endif
