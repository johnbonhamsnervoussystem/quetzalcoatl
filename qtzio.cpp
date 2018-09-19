/* Routines for Quetz I/O */
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
#include "qtzio.h"
#include "tei.h"
#include "time_dbg.h"
using namespace Eigen ;
using namespace std ;

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
  std::ifstream jobfile ;
  std::vector<double> atnum ;
  std::vector<std::vector<double>> t_c ;
  time_dbg read_input_time = time_dbg("read_input") ;

  jobfile.open( inpfile, std::ifstream::in ) ;

  while( getline( jobfile, line)){
    if ( line.substr(0,7) == "jobtyp " ){
      i_junk = stoi(line.substr(9)) ;
    }  else if ( line.substr(0,7) == "basis  " ){
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

void getmel( string file1, string file2, vector<tei>& intarr, common& com) {
/* 
 * Read matrix elements from phfinp .
 *
 * */
  int nbasis=0 ;
  int i ;
  int j ;
  int i_v ;
  int j_v ;
  int k_v ;
  int l_v ;
  double int_v ;
  MatrixXd ovl ;
  MatrixXd ham ;
  ifstream input_file ;
  string line ;
  tei tmp_tei ;
  
  /* This reads in data formatted from a matrix element file 
   to a second file and finally to here. */
  nbasis = com.nbas() ;

  ham.resize( nbasis, nbasis) ;
  ovl.resize( nbasis, nbasis) ;
  input_file.open(file1) ;
  
  if ( input_file.is_open() ) {

    getline( input_file, line) ;
    getline( input_file, line) ;
  
    if( line != " OVERLAP") {
      cout << "No overlap in the input file." << endl ;
      }
  
    for( i = 0 ; i < nbasis ; i++ ){
      for( j = 0 ; j <= i ; j++ ){
        getline( input_file, line) ;
        ovl(i,j) = stof(line) ;
        ovl(j,i) = stof(line) ;
        }
      }
  
    getline( input_file, line) ;
  
    if( line != " CORE HAMILTONIAN") {
      cout << "No core hamiltonian." << endl ;
      }
  
    for( i = 0 ; i < nbasis ; i++ ){
      for( j = 0 ; j <= i ; j++ ){
        getline( input_file, line) ;
        ham(i,j) = stof(line) ;
        ham(j,i) = stof(line) ;
        }
      }
  
    input_file.close() ;
    com.setS( ovl) ;
    ovl.resize(0,0) ;
    com.setH( ham) ;
    ham.resize(0,0) ;

    }

    input_file.open(file2) ;

    if ( input_file.is_open() ) {
      while( getline( input_file, line) ){ 
        i_v = stoi(line.substr(1,2)) ;
        j_v = stoi(line.substr(4,5)) ;
        k_v = stoi(line.substr(7,8)) ;
        l_v = stoi(line.substr(10,11)) ;
        int_v = stof(line.substr(13,25)) ;
        tmp_tei.set( i_v, j_v, k_v, l_v, int_v) ;
        intarr.push_back(tmp_tei) ;
        }
      }

  input_file.close() ;


  return ;

} 

void rdsdet ( int nbasis, vector<string>& matel, vector<hfwfn>& det) {
/* 
 * This routine gets wavefunctions.
 *
 * This routine does no dimensioning for you.  Pass it preallocated stuff please.
 *
 * */
  int c_len ;
  int nfiles=matel.size() ;
  int jcol ;
  int jrow ;
  ifstream i_file ;
  string file_read ;
  string line ;
  string rc_str ;
  string ic_str ;
  double rc_fl ;
  double ic_fl ;
  Eigen::MatrixXcd moc;

  c_len = nbasis*nbasis*4 ;
  moc.resize( 2*nbasis, 2*nbasis) ;

  for ( int i=0; i < nfiles; i++ ){

    file_read = matel[i] ;
    i_file.open(file_read) ;

    if ( i_file.is_open() ) {

      for ( int j=0; j < c_len; j++ ){

        jrow = j % (2*nbasis) ;
        jcol = j/(2*nbasis) ;
        getline( i_file, line) ;
        rc_str = line.substr(0,14) ;
        ic_str = line.substr(16,30) ;
        rc_fl = stof(rc_str) ;
        ic_fl = stof(ic_str) ;
        moc( jrow, jcol) = cd (rc_fl, ic_fl) ;
      }

    }

    det[i].set_mos( moc) ;
    i_file.close() ;

  }

  return ;

}

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

template void read_eigen_bin(const Eigen::VectorXf&, std::ofstream&) ;
template void read_eigen_bin(const Eigen::VectorXd&, std::ofstream&) ;
template void read_eigen_bin(const Eigen::VectorXcf&, std::ofstream&) ;
template void read_eigen_bin(const Eigen::VectorXcd&, std::ofstream&) ;
template void read_eigen_bin(const Eigen::MatrixXf&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXcf&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXcd&, std::ifstream&) ;

