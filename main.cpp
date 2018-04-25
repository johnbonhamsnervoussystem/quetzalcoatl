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

int main() {
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
  int nbasis ;
  int nbas ;
  int nele ;
  int nalp ;
  int nbet ;
  float fjunk ;
  hfwfn detrr ;
  hfwfn detcr ;
  hfwfn detru ;
  hfwfn detcu ;
  hfwfn detrg ;
  hfwfn detcg ;

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
    
  nbasis = com.nbas() ;
  nele = com.nele() ;
  nalp = com.nalp() ;
  nbet = com.nbet() ;

  getmel( "f00.fi1s", "f00.fi2s", intarr, com) ;

//  std::cout << "real restricted" << std::endl ;
//  detrr.set( com, intarr, "rrhf") ;
//  detrr.prt_mos() ;
//  detrr.prt_ene( com.nrep()) ;
//  std::cout << std::endl << "real unrestricted" << std::endl ;
//  detru.set( com, intarr, "ruhf") ;
//  detru.prt_mos() ;
//  detru.prt_ene( com.nrep()) ;
  std::cout << std::endl << "real general" << std::endl ;
  detrg.set( com, intarr, "rghf") ;
  detrg.prt_mos() ;
  detrg.prt_ene( com.nrep()) ;
//  std::cout << std::endl << "cmplx restricted" << std::endl ;
//  detcr.set( com, intarr, "crhf") ;
//  detcr.prt_mos() ;
//  detcr.prt_eig() ;
//  detcr.prt_ene( com.nrep()) ;
//  std::cout << std::endl << std::endl ;
//  std::cout << std::endl << "cmplx unrestricted" << std::endl ;
//  detcu.set( com, intarr, "cuhf") ;
//  detcu.prt_mos() ;
//  detcu.prt_eig() ;
//  detcu.prt_ene( com.nrep()) ;
  std::cout << std::endl << "cmplx general" << std::endl ;
  detcg.set( com, intarr, "cghf") ;
  detcg.prt_mos() ;
  detcg.prt_eig() ;
  detcg.prt_ene( com.nrep()) ;
  
  return 0 ;

} /* End Quetzacoatl */
