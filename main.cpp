/* Primary routine for Quetzalcoatl */
#include<iostream>
#include<complex>
#include<string>
#include<fstream>
#include<vector>
#include<Eigen/Dense>
#include<fstream>
#include "tei.h"
#include "qtzio.h"
#include "common.h"
#include "hfwfn.h"
#include "util.h"
#include "evalm.h"

int main() {
/*
 * To Do :: 
 *  
 *   implement stability 
 *   matrix elements between determinants
 *   wigner matrix
 *   Eventually, I will need to write a routine which grabs a chunk of memory,
 *   rather than constantly allocating and dealllocating.  For now, I need things
 *   that work and produce useful results.  
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
  hfwfn detrr1 ;
  Eigen::MatrixXf s ;
  Eigen::MatrixXf h ;
  Eigen::MatrixXf h_full ;
  Eigen::MatrixXf s_full ;
  Eigen::MatrixXf mos ;
  Eigen::MatrixXf mosf ;
  Eigen::MatrixXf tmp ;

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

  getmel( "f00.fi1s", "f00.fi2s", intarr, com) ;

  h_full.resize( 2*com.nbas(), 2*com.nbas()) ;
  s_full.resize( 2*com.nbas(), 2*com.nbas()) ;
  h.resize( com.nbas(), com.nbas()) ;
  s.resize( com.nbas(), com.nbas()) ;
  mos.resize( com.nbas(), 2*com.nbas()) ;
  mosf.resize( 2*com.nbas(), 2*com.nbas()) ;
  tmp.resize( 2*com.nbas(), 2*com.nbas()) ;
  h = com.getH() ;
  s = com.getS() ;
  oao ( com.nbas(), h, s) ;
  s_full.setZero() ;
  s_full.block( 0, 0, com.nbas(), com.nbas()) = s ;
  s_full.block( com.nbas(), com.nbas(), com.nbas(), com.nbas()) = s ;

  detrr1.init( com, intarr, "ruhf") ;
  detrr1.prt_ene ( com.nrep()) ;
  std::cout << " Retrieving the mos" << std::endl;
  detrr1.get_mos ( mos) ;
  mosf.setZero() ;
  mosf.block( 0, 0, com.nbas(), com.nbas()) = mos.block( 0, 0, com.nbas(), com.nbas()) ;
  mosf.block( com.nbas(), com.nbas(), com.nbas(), com.nbas()) = mos.block( 0, com.nbas(), com.nbas(), com.nbas()) ;
  std::cout << " CtSC " << std::endl;
  tmp = mosf.adjoint()*s_full*mosf ;
  std::cout << tmp << std::endl ;

  oao( com.nbas(), detrr1, s) ;
  detrr1.get_mos ( mos) ;
  mosf.setZero() ;
  mosf.block( 0, 0, com.nbas(), com.nbas()) = mos.block( 0, 0, com.nbas(), com.nbas()) ;
  mosf.block( com.nbas(), com.nbas(), com.nbas(), com.nbas()) = mos.block( 0, com.nbas(), com.nbas(), com.nbas()) ;
  std::cout << " Ct*C " << std::endl;
  tmp = mosf.adjoint()*mosf ;
  std::cout << tmp << std::endl ;

  h_full.setZero() ;
  h_full.block( 0, 0, com.nbas(), com.nbas()) = h ;
  h_full.block( com.nbas(), com.nbas(), com.nbas(), com.nbas()) = h ;

  oao( com.nbas(), intarr, tmparr, s) ;

  fjunk = fockop( com, h_full, tmparr, detrr1, detrr1) ;
  std::cout << fjunk + com.nrep() << std::endl ;

  return 0 ;

} /* End Quetzacoatl */

