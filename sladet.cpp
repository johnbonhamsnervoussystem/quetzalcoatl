#include <complex>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "qtzio.h"
#include "sladet.h"
#include <sys/stat.h>

/**
 *
 Check if a file exists
 * @return true if and only if the file exists, false else
 */

const std::string wfnIO = "qtz.wfn.bin" ;

template<typename s, int r, int c>
void save_slater_det( sladet< s, r, c>& w, int cntl){
/* A controlling routine to know where to create a new external file
   or to append to an existing file.

  cntl : 
    0 - (default) write to a new file or overwrite one that already exists 
    1 - Append to an existing file.  If it doesn't exist, create it.       */
  bool ioerr = false ;
  std::ofstream F_OUT ;
  struct stat buf;

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
    std::cout << "Unrecognized option in save_slater_det" << std::endl ;
    std::cout << "Exiting Quetzalcoatl" << std::endl ;
    exit(EXIT_FAILURE) ;
  }

  if ( ioerr ) {
    std::cout << "IOERR" << std::endl ;
    std::cout << "Exiting Quetzalcoatl" << std::endl ;
    exit(EXIT_FAILURE) ;
    }

  /* Presumably we have a file to write to now! */
  write_to_bin( w, F_OUT) ;

  F_OUT.close() ;

  return ;

}

template void save_slater_det( sladet< double, Eigen::Dynamic, Eigen::Dynamic>&, int) ;
template void save_slater_det( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&, int) ;

template<typename s, int r, int c>
void write_to_bin( sladet< s, r, c>& w, std::ofstream& F_OUT){ 
/* MO coefficients and eigenvalues to an external file. */

  F_OUT << w.e_scf ;
  write_eigen_bin( w.moc, F_OUT) ; 
  write_eigen_bin( w.eig, F_OUT) ;

  F_OUT.close() ;

  return ;
  
} ; 

template void write_to_bin( sladet< double, Eigen::Dynamic, Eigen::Dynamic>&, std::ofstream&) ;
template void write_to_bin( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&, std::ofstream&) ;

template<typename s, int r, int c>
void load_slater_det( sladet< s, r, c>& w, int cntl) {
/* A controlling routine to read slater determinants from file.  It's super ugly and needs tlc at some point

  cntl : 
    0 - (default) read the first SlaDet. 
    N - Read N Slater Determinants 
  */
  bool ioerr = false ;
  std::ifstream F_IN ;
  struct stat buf;

    /* Check if the file exists */
  if ( stat(wfnIO.c_str(), &buf) == 0) {
      try {
        F_IN.open( wfnIO, std::ifstream::binary | std::ifstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
  } else {
    std::cout << "File does not exist." << std::endl ;
    std::cout << "Exiting Quetzalcoatl" << std::endl ;
    exit(EXIT_FAILURE) ;
  }

  if ( ioerr ) {
    std::cout << "IOERR" << std::endl ;
    std::cout << "Exiting Quetzalcoatl" << std::endl ;
    exit(EXIT_FAILURE) ;
    }

  /* Presumably we have a file to write to now! */
  read_sladet( w, F_IN) ;

  F_IN.close() ;

  return ;

} ;

template void load_slater_det( sladet< double, Eigen::Dynamic, Eigen::Dynamic>&, int) ;
template void load_slater_det( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&, int) ;

template<typename s, int r, int c>
void read_sladet( sladet< s, r, c>& w, std::ifstream& F_IN){ 
/* Write the wfntype, MO coefficients and eigenvalues to an external file. */

  F_IN >> w.e_scf ;
  read_eigen_bin( w.moc, F_IN) ;
  read_eigen_bin( w.eig, F_IN) ;

  F_IN.close() ;

  return ;

} ;

template void read_sladet( sladet< double, Eigen::Dynamic, Eigen::Dynamic>&, std::ifstream&) ;
template void read_sladet( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&, std::ifstream&) ;
