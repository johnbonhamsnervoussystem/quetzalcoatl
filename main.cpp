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

int main() {
/*
 * To Do :: 
 *  
 *   implement stability 
 *   matrix elements between determinants
 *   wigner matrix
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

  typedef std::complex<float> cf ;
  int junk ;
  int job=0  ;
  float fjunk ;
  hfwfn detrr ;
  hfwfn detru ;
  hfwfn detrg ;
  Eigen::MatrixXf s ;

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

  detrr.init( com, intarr, "rrhf" ) ;
  detru.init( com, intarr, "ruhf" ) ;
  detrg.init( com, intarr, "rghf" ) ;
  s.resize( com.nbas(), com.nbas()) ;
  s = com.getS() ;
  oao ( com.nbas(), detrr, s) ;
  oao ( com.nbas(), detru, s) ;
  oao ( com.nbas(), detrg, s) ;
  s.resize( 0, 0) ;
  
  return 0 ;

} /* End Quetzacoatl */
