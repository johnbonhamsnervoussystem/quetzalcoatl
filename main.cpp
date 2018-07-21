/* Primary routine for Quetzalcoatl */
#include "constants.h"
#include <complex>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <sstream>
#include <fstream>
#include <ctime>
#include <libint2.h>
#include "binio.h"
#include "common.h"
#include "evalm.h"
#include "hfwfn.h"
#include "integr.h"
#include "project.h"
#include "qtzio.h"
#include "solver.h"
#include "tei.h"
#include "util.h"
#include "wigner.h"

int main(int argc, char *argv[]) {
/*
 * This is my own implementation of various electronic stucture methods.  
 * There is lots of work to do and improvements to be made and the purpose
 * is supposed to be minamly pedagoical as well as producing results. 
 *
 * To Do :: 
 *   - Parse options and set route through program.
 *   -implement stability 
 *   -matrix elements between determinants
 *   -Eventually, I will need to write a routine which grabs a chunk of memory,
 *   rather than constantly allocating and dealllocating.  For now, I need things
 *   that work and produce useful results.  
 *   - Debug spin rotations on a determinant.
 *   - Interface with non-gaussian routines for matrix elements
 *
 * */

/* 
  job_file - file to read information about the job
  readfile - name of the job file being read in 
  line - holds information being read in
  job_dim - A class to carry relevant dimensioning information through
	the program.
  
  ovl - matrix containing the overlp matrix
  ham - matrix containing the core hamiltonian
  xmat - matrix containing the transformation to an orthogonal ao basis.
  */
  common com ;
  std::vector<tei> intarr ;
  std::vector<tei> tmparr ;

  int job=0  ;
  double cx = 0.0e0 ;
  double cy = 0.0e0 ;
  double cz = 0.0e0 ;
  double r = 0.0e0 ;
  double r2 = 0.0e0 ;
  double n_rep = 0.0e0 ;
  cd ejunk ;
  cd ojunk ;
  hfwfn det1 ;
  hfwfn det2 ;
  std::vector<hfwfn> tst_vec ;
  std::vector<std::string> wfn_vec ;
  std::string trden="tden.rwf" ;
  std::string fokmat="fmat.rwf" ;
  Eigen::MatrixXd c ;
  Eigen::VectorXi a ;
  /* File reading and header variables. */
  std::stringstream ss ;
  std::string inpfile ;
  std::time_t t = std::time(0) ;
  std::tm* now = std::localtime(&t) ;
/*Eigen::MatrixXd mosf ;
  Eigen::MatrixXd tmp ; */

/* Open and read the input file */

/*   Let's just us a test file for now.
 *   std::cin >> readfile ; */
  ss << argv[1] ;
  ss >> inpfile ;
  std::cout << " ---   ---   ---   --- " << std::endl << "|" << std::endl 
    << "|  Quetzalcoatl v1.0  " << std::endl ;
  std::cout << "| " << std::endl ;
  std::cout << " ---   ---   ---   --- " << std::endl << "|" << std::endl ;
  std::cout << "|  Started at  " << std::endl << "| " << now->tm_hour  << " : " 
    << now->tm_min + 1 << " : " << now->tm_sec << std::endl ;
  std::cout << "|  on  " << std::endl << "| " << now->tm_year + 1900 << " - " 
    << now->tm_mon + 1 << " - " << now->tm_mday << std::endl ;
  std::cout << "| " << std::endl ;
  std::cout << " ---   ---   ---   --- " << std::endl << "|" << std::endl ;
  std::cout << "|  Reading from input : " << inpfile << std::endl ;
  std::cout << "|  Basis Set : " << std::cout << com.bnam()  << std::endl ;
  std::cout << "| " << std::endl ;

  read_input( com, inpfile) ;
  std::cout << "| " << std::endl ;
  std::cout << " ---   ---   ---   --- " << std::endl ;

  c.resize( com.natm(), 3) ;
  a.resize( com.natm()) ;
  c = com.getC() ;
  a = com.getA() ;
  

  for ( int i = 0; i < com.natm(); i++) {
    for ( int j = i+1; j < com.natm(); j++) {
      cx = c( i, 0) - c( j, 0) ;
      cy = c( i, 1) - c( j, 1) ;
      cz = c( i, 2) - c( j, 2) ;
      r2 = cx*cx + cy*cy + cz*cz ;
      r = sqrt(r2) ;
      n_rep += static_cast<double>(a(i))*static_cast<double>(a(j))/r ;
      }
    }

  com.nrep( n_rep ) ;

  LIBINT2_PREFIXED_NAME(libint2_static_init)() ;
  LIBINT2_PREFIXED_NAME(libint2_static_cleanup)() ;
//  getmel( "./test_fil/f00.fi1s", "./test_fil/f00.fi2s", intarr, com) ;
//  s.resize( com.nbas(), com.nbas()) ;
//  s = com.getS() ;
//
//  h.resize( com.nbas(), com.nbas()) ;
//  h = com.getH() ;
//  oao( com.nbas(), h, s) ;
//  h_full.resize( 2*com.nbas(), 2*com.nbas()) ;
//  h_full.setZero() ;
//
//  det1.fil_mos( com.nbas(), h_full, 6) ;
//  det2.fil_mos( com.nbas(), h_full, 6) ;
//
//  h_full.block( 0, 0, com.nbas(), com.nbas()).real() = h ;
//  h_full.block( com.nbas(), com.nbas(), com.nbas(), com.nbas()).real() = h ;
//
//  wfn_vec.push_back("./test_fil/wfn1.det") ;
//  wfn_vec.push_back("./test_fil/wfn2.det") ;
//
//  tst_vec.push_back( det1) ;
//  tst_vec.push_back( det2) ;
//
//  rdsdet( com.nbas(), wfn_vec, tst_vec) ;
//
//  oao( com.nbas(), tst_vec[0], s) ;
//  oao( com.nbas(), tst_vec[1], s) ;
//  oao( com.nbas(), intarr, tmparr, s) ;
//
//  trci( com, tst_vec, h_full, tmparr, trden, fokmat) ;
//
  return 0 ;

} /* End Quetzacoatl */

