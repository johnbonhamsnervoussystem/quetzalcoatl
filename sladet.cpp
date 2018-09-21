#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "qtzio.h"
#include "sladet.h"
#include <sys/stat.h>

/**
 * Check if a file exists
 * @return true if and only if the file exists, false else
 */
bool fileExists(const std::string& file) {


const std::string wfnIO = "qtz.wfn.bin" ;

template<typename s, int r, int c>
void save_slater_det( sladet< s, r, c>& w, int& cntl){
/* A controlling routine to know where to create a new external file
   or to append to an existing file.

  cntl : 
    0(default) - write to a new file or overwrite one that already exists */
  std::ofstream F_OUT ;
  struct stat buf;

  /* Check if the file exists. */
  if ( stat(wfnIO.c_str(), &buf) == 0) {

    try {
      F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
    } catch ( std::fstream::failure& e) {
      std::cout << e.what() << std::endl ;
    } catch (...) {
      std::cout << " Default excepetion inside write_to_bin " << std::endl ;
      }

  return ;

}

template<typename s, int r, int c>
void write_to_bin( sladet< s, r, c>& w, std::ofstream& F_OUT){ 
/* Write the wfntype, MO coefficients and eigenvalues to an external file. */

  F_OUT.write( s.wfntp) ; F_OUT.write( s.e_scf) ;
  write_eigen_bin( s.moc, F_OUT) ; write_eigen_bin( s.eig, F_OUT) ;

  F_OUT.close() ;

  return ;
  
} ; 

template void write_to_bin( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&) ;

template<typename s, int r, int c>
void load_slater_det( sladet< s, r, c>& w) {
/* */
  std::ofstream F_IN ;

  try {
    F_IN.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
  } catch ( std::fstream::failure& e) {
    std::cout << e.what() << std::endl ;
  } catch (...) {
    std::cout << " Default excepetion inside load_sd " << std::endl ;
    }

  F_IN.read( s.e_scf) ;
  read_eigen_bin( s.moc, F_IN) ; read_eigen_bin( s.eig, F_IN) ;

  F_IN.close() ;

  return ;

} ;

template void load_sd( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&) ;

