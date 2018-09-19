#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "qtzio.h"
#include "sladet.h"

const std::string wfnIO = "qtz.wfn.bin" ;

template<typename s, int r, int c>
void write_to_bin( sladet< s, r, c>& w){ 
/* Write the wfntype, MO coefficients and eigenvalues to an external file. */
  std::ofstream F_OUT ;
  
  try {
    F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
  } catch ( std::fstream::failure& e) {
    std::cout << e.what() << std::endl ;
  } catch (...) {
    std::cout << " Default excpetion inside write_to_bin " << std::endl ;
    }

  F_OUT.write( s.wfntp) ; F_OUT.write( s.e_scf) ;
  write_eigen_bin( s.moc, F_OUT) ; write_eigen_bin( s.eig, F_OUT) ;

  F_OUT.close() ;

  return ;
  
} ; 

template void write_to_bin( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&) ;

template<typename s, int r, int c>
void load_sd( sladet< s, r, c>& w) {
/* Write the wfntype, MO coefficients and eigenvalues to an external file. */
  std::ofstream F_IN ;

  try {
    F_IN.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
  } catch ( std::fstream::failure& e) {
    std::cout << e.what() << std::endl ;
  } catch (...) {
    std::cout << " Default excpetion inside load_sd " << std::endl ;
    }

  F_IN.read( s.e_scf) ;
  read_eigen_bin( s.moc, F_IN) ; read_eigen_bin( s.eig, F_IN) ;

  F_IN.close() ;

  return ;

} ;

template void load_sd( sladet< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&) ;

