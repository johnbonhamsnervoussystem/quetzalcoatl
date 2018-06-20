/* Primary routine for Quetzalcoatl */
#include "constants.h"
#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
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

int main() {
/*
 * This is my own implementation of various electronic stucture methods.  
 * There is lots of work to do and improvements to be made and the purpose
 * is supposed to be minaly pedagoical as well as producing results.  More 
 * accurate and faster implementations will always be possible.  This is 
 * more like a sandbox to try new things.
 *
 * To Do :: 
 *  
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
   
  std::ifstream job_file ;
  std::string readfile ;
  std::string line ;
  common com ;
  std::vector<tei> intarr ;
  std::vector<tei> tmparr ;

  typedef std::complex<float> cf ;
  int junk ;
  int job=0  ; 
  float fjunk ;
  cf cjunk ;
  cf ojunk ;
  hfwfn det1 ;
  std::vector<hfwfn> tst_vec ;
  Eigen::MatrixXf h ;
  Eigen::MatrixXf s ;
  Eigen::MatrixXcf h_full ;
  Eigen::MatrixXf s_full ;
  Eigen::MatrixXcf mos ;
/*Eigen::MatrixXf mosf ;
  Eigen::MatrixXf tmp ; */

/* Open and read the input file */

/*   Let's just us a test file for now.
 *   std::cin >> readfile ; */
  readfile = "test_file.qtz" ;
  job_file.open(readfile) ;
  while( getline( job_file, line)){
    if ( line.substr(0,7) == "jobtyp " ){
      job = stoi(line.substr(9)) ;
    }  else if ( line.substr(0,7) == "nbasis " ){
      junk = stoi(line.substr(9)) ;
      com.nbas(junk) ;
    }  else if ( line.substr(0,7) == "nn rep " ) {
      fjunk = stof(line.substr(9)) ;
      com.nrep( fjunk) ;
    }  else if ( line.substr(0,7) == "n2ei   " ) {
      junk = stoi(line.substr(9)) ;
      com.ntei(junk) ;
    }  else if ( line.substr(0,7) == "natoms " ) {
      junk = stoi(line.substr(9)) ;
      com.natm(junk) ;
    }  else if ( line.substr(0,7) == "nelec  " ) {
      junk = stoi(line.substr(9)) ;
      com.nele(junk) ;
    }  else if ( line.substr(0,7) == "nalpha " ) {
      junk = stoi(line.substr(9)) ;
      com.nalp(junk) ;
    }  else if ( line.substr(0,7) == "nbeta  " ) {
      junk = stoi(line.substr(9)) ;
      com.nbet(junk) ;
    }
  }

  readcx( com.nbas()) ;
/*  getmel( "f00.fi1s", "f00.fi2s", intarr, com) ;
  det1.init( com, intarr, "cghf") ;

  s.resize( com.nbas(), com.nbas()) ;
  s = com.getS() ;

  h.resize( com.nbas(), com.nbas()) ;
  h = com.getH() ;

  oao( com.nbas(), det1, s) ;
  oao( com.nbas(), h, s) ;
  oao( com.nbas(), intarr, tmparr, s) ;

  h_full.resize( 2*com.nbas(), 2*com.nbas()) ;
  h_full.setZero() ;
  h_full.block( 0, 0, com.nbas(), com.nbas()).real() = h ;
  h_full.block( com.nbas(), com.nbas(), com.nbas(), com.nbas()).real() = h ;

  tst_vec.push_back(det1) ;
  trci( com, tst_vec, h_full, tmparr) ; */

  return 0 ;

} /* End Quetzacoatl */

