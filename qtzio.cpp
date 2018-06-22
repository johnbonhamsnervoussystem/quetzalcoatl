/* Routines for Quetz I/O */
#include "constants.h"
#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <sys/stat.h>
#include "qtzio.h"
#include "tei.h"
#include "common.h"
using namespace Eigen ;
using namespace std ;


void getmel( string file1, string file2, vector<tei>& intarr, common& com) {
/* 
 * Read matrix elements from phfinp .
 *
 * */
  int nbasis=0 ;
  int nbsuse=0 ;
  int i ;
  int j ;
  int i_v ;
  int j_v ;
  int k_v ;
  int l_v ;
  float int_v ;
  MatrixXf ovl ;
  MatrixXf ham ;
  ifstream input_file ;
  string line ;
  tei tmp_tei ;
  
  /* This reads in data formatted from a matrix element file 
   to a second file and finally to here. */
  nbasis = com.nbas() ;
  nbsuse = com.nbsu() ;


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
    com.setS( ovl, nbasis) ;
    ovl.resize(0,0) ;
    com.setH( ham, nbasis) ;
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
  float rc_fl ;
  float ic_fl ;
  Eigen::MatrixXcf moc;

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
        moc( jrow, jcol) = cf (rc_fl, ic_fl) ;
      }

    }

    det[i].set_mos( moc) ;
    i_file.close() ;

  }

  return ;

}

void w_eigen_binary( ) {
/* Passed an eigen matrix and ofstream, write the contents. */

}

