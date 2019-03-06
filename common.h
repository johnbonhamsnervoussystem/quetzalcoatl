#include <Eigen/Dense>
#include "tei.h"
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
    hub_opt - Options for hubbard

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
    U_interact - interaction term. U for Hubbard, G for Pairing

  Algorithm Control :
    max_scf_iter - the maximum number of scf iterations
    max_pn_iter - the maximum number of particle number iterations for HFB
    scf_convergence_threshold - the threshold at which we consider the wavefunction
      converged.
    print - whether to do debug print or not
    level_shift - Whether to use level shifting
    lshift - an initial value for level shifting

  Stored Matrices
   s - overlap matrix
   h - core hamiltonian
   a - atomic number
   coord - coordinates matrix

  Projection Data
    nbr_g - particle number projection grid
*/

private :

  int hamiltonian ;
  int method ;
  int hubbard ;
  int *hub_dim ;
  int nbasis ;
  int n2ei ;
  int natoms ;
  int nel ;
  int nal ;
  int nbe ;
  double chemical_potential ;
  double U_interact ;
  double nn ;
  std::string basis_name ;

  int max_scf_iter ;
  int max_pn_iter ;
  int print ;
  bool level_shift ;
  double lshift ;
  double scf_convergence_threshold ;

  Eigen::MatrixXd s_c ;
  Eigen::MatrixXd xs_c ;
  Eigen::MatrixXd h_c ;
  Eigen::VectorXd a_c ;
  Eigen::MatrixXd coord ;

  std::vector<tei> r12int ;
  int nbr_g ;

public :
/* Initializers to set default values */
  common( void) ;

/* Set the data */
  void hamil( int n) ;
  void methd( int n) ;
  void hub_opt( int n) ;
  void hub_n( int n, int i) ;
  void nbas( int n) ;
  void ntei( int n) ;
  void natm( int n) ;
  void nele( int n) ;
  void nalp( int n) ;
  void nbet( int n) ;
  void mu( double n) ;
  void nrep( double f) ;
  void bnam( std::string n) ;
  void ngrid( int n) ;

/* Algorithm control*/
  void mxscfit( int n) ;
  void mxpnit( int n) ;
  void prt( int n) ;
  void scfthresh( double d) ;
  void lvlshft( double d) ;
  bool use_shift( void) ;

/* Coordinates */
  void setcoord( std::vector<std::vector<double>> c) ;

/* Set matrix elements */
  void setS( Eigen::MatrixXd s) ;
  void setXS( Eigen::MatrixXd xs) ;
  void setH( Eigen::MatrixXd h) ;
  void setA( std::vector<double> a) ;
  void setC( std::vector<std::vector<double>> c) ;
  void setr12( std::vector<tei> intarr) ;
  void setU( double u) ;

/* Retrieve the data */
  int hamil( void) ;
  int methd( void) ;
  int hub_opt( void) ;
  int hub_n( int i) ;
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
  int mxpnit( void) ;
  int prt( void) ;
  double scfthresh( void) ;
  double lvlshft( void) ;

/* Retrieve a matrix */
  Eigen::MatrixXd getS( void) ;
  Eigen::MatrixXd getXS( void) ;
  Eigen::MatrixXd getH( void) ;
  Eigen::VectorXd getA( void) ;
  Eigen::MatrixXd getC( void) ;
  void getr12( std::vector<tei>*& intarr) ;
  double getU( void) ;


} ;
#endif
