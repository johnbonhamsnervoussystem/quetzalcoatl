/* Routines for Quetz I/O */
#include <algorithm>
#include "common.h"
#include <complex>
#include "constants.h"
#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <sys/stat.h>
#include "qtzcntrl.h"
#include "qtzio.h"
#include "tei.h"
#include "time_dbg.h"

void read_input( common& com, const std::string& inpfile){
 /* Read the input file and save the information into common. */
  int i_junk ;
  int at_count = 0 ;
  std::vector<double> tmp(3) ;
  double d_junk ;
  size_t pos = 0 ;
  std::string line ;
  std::string s_junk ;
  std::string delim = "," ;
  std::string::iterator end_pos ;
  std::ifstream jobfile ;
  std::vector<double> atnum ;
  std::vector<std::vector<double>> t_c ;
  time_dbg read_input_time = time_dbg("read_input") ;

  jobfile.open( inpfile, std::ifstream::in ) ;

  while( getline( jobfile, line)){
    if ( line.substr(0,7) == "hamilt " ){
      s_junk = line.substr(9) ;
      strip_lower( s_junk) ;
/*
      Molecular/mol
      Semi-Empirical/se
*/
      if ( s_junk == "mol"){
/* Molecular Electronic Hamiltonian */
        com.hamil(1) ;
      } else if ( s_junk == "se"){
/* Semi-Empirical Hamiltonian */
        std::cout << "  NOT IMPLEMENTED " << std::endl ;
        com.hamil(1) ;
      } else {
        com.hamil(1) ;
        }

    } else if ( line.substr(0,7) == "method " ){
/*
        real/r
        complex/c

        restricted/r
        unrestricted/u
        genralized/g

	Hartree-Fock/HF
        Hartree-Fock-Bogliubov/HFB
*/
      s_junk = line.substr(9) ;
      strip_lower( s_junk) ;
      if ( s_junk.substr(0,2) == "rr" ) { com.methd(1) ;}
      else if ( s_junk.substr(0,2) == "cr" ) { com.methd(2) ;}
      else if ( s_junk.substr(0,2) == "ru" ) { com.methd(3) ;}
      else if ( s_junk.substr(0,2) == "cu" ) { com.methd(4) ;}
      else if ( s_junk.substr(0,2) == "rg" ) { com.methd(5) ;}
      else if ( s_junk.substr(0,2) == "cg" ) { com.methd(6) ;}
      if ( s_junk.substr(2) == "hf" ) { com.methd(10) ;}
      else if ( s_junk.substr(2) == "hfb" ) { com.methd(20) ;}
    } else if ( line.substr(0,7) == "basis  " ){
      s_junk = line.substr(9) ;
      com.bnam( s_junk) ;
    }  else if ( line.substr(0,7) == "geom   " ) {
      getline( jobfile, line) ;
      while( line.substr(0,5) != " end" ){
        at_count++ ;
        pos = line.find(delim) ;
        d_junk = stod(line.substr( 0, pos)) ;
        atnum.push_back( d_junk) ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[0] = stod(line.substr( 0, pos)) ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[1] = stod(line.substr( 0, pos)) ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[2] = stod(line.substr( 0, pos)) ;
        line.erase( 0, pos + delim.length()) ;
        getline( jobfile, line) ;
        t_c.push_back( tmp) ;
        }
      com.natm( at_count) ;
      com.setA( atnum) ;
      com.setC( t_c) ;
    }  else if ( line.substr(0,7) == "nalpha " ) {
      i_junk = stoi(line.substr(9)) ;
      com.nalp(i_junk) ;
    }  else if ( line.substr(0,7) == "nbeta  " ) {
      i_junk = stoi(line.substr(9)) ;
      com.nbet(i_junk) ;
    }
  }

  com.nele(com.nalp() + com.nbet()) ;
  read_input_time.end() ;

  return ;

} 

bool open_binary( std::ofstream& F_OUT, int cntl) {
/*
  cntl :
    0 - (default) write to a new file or overwrite one that already exists
    1 - Check if a file exists.  If not create it.
*/
  const std::string wfnIO = "qtz.wfn.bin" ;
  bool ioerr = false ;
  struct stat buf ;

  if ( cntl == 0 ){
    try {
      F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
    } catch ( std::fstream::failure& e) {
      std::cout << e.what() << std::endl ;
      ioerr = true ;
    } catch (...) {
      std::cout << " Default excepetion inside write_to_bin " << std::endl ;
      ioerr = true ;
      }
  } else if ( cntl == 1 ){
    /* Check if the file exists */
    if ( stat(wfnIO.c_str(), &buf) == 0) {
      try {
        F_OUT.open( wfnIO, std::ofstream::app | std::ofstream::binary | std::ofstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
    } else {
      try {
        F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
    }
  } else {
    qtzcntrl::shutdown( "Unrecognized Option (open_binary)") ;
    }

  return ioerr ;

} ;

bool open_binary( std::ifstream& F_IN) {
/*
  (default) Check if a file exists and if it does open it
*/
  const std::string wfnIO = "qtz.wfn.bin" ;
  bool ioerr = false ;
  struct stat buf ;

  /* Check if the file exists */
  if ( stat(wfnIO.c_str(), &buf) == 0) {
    try {
      F_IN.open( wfnIO, std::ofstream::binary | std::ofstream::in) ;
    } catch ( std::fstream::failure& e) {
      std::cout << e.what() << std::endl ;
      ioerr = true ;
    } catch (...) {
      std::cout << " Default exception inside write_to_bin " << std::endl ;
      ioerr = true ;
      }
  } else {
    qtzcntrl::shutdown( "FIle does not exist (open_binary)") ;
    }

  return ioerr ;

} ;

void strip_lower( std::string& s) {
/*
  Wrap together whitespace trimming and lower casing for input file parsing
*/
  std::string::iterator end_pos ;
  end_pos = std::remove( s.begin(), s.end(), ' ') ;
  s.erase( end_pos, s.end()) ;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower) ;

  return ;

} ;

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ofstream& F_OUT) {
  typename matrix::Index rows=m.rows(), cols=m.cols();

  F_OUT.write( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;

  return ;

}

template void write_eigen_bin(const Eigen::VectorXf&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::VectorXd&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::VectorXcf&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::VectorXcd&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::MatrixXf&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::MatrixXd&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::MatrixXcf&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::MatrixXcd&, std::ofstream&) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::ifstream& F_IN) {
  typename matrix::Index rows=m.rows(), cols=m.cols();

  F_IN.read( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;

  return ;

}

template void read_eigen_bin(const Eigen::VectorXf&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::VectorXd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::VectorXcf&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::VectorXcd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXf&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXcf&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXcd&, std::ifstream&) ;

