#include <complex>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "qtzcntrl.h"
#include "qtzio.h"
#include "wfn.h"
#include <sys/stat.h>

const std::string wfnIO = "qtz.wfn.bin" ;

template<typename s, int r, int c>
void save_wfn( wfn< s, r, c>& w, int cntl){
/*
  Save a slater determinant to a binary file for passing between
  various routines.

  cntl :
    0 - (default) write to a new file or overwrite one that already exists
    1 - Append to an existing file.  If it doesn't exist, create it.
*/
  std::ofstream F_OUT ;

  if ( open_binary( F_OUT, cntl) ) {
    qtzcntrl::shutdown( "IO Error (save_wfn)") ;
    }

  F_OUT << w.wfntyp ;
  F_OUT << w.e_scf ;
  write_eigen_bin( w.moc, F_OUT) ;
  write_eigen_bin( w.eig, F_OUT) ;

  F_OUT.close() ;

  return ;

}

template void save_wfn( wfn< double, Eigen::Dynamic, Eigen::Dynamic>&, int) ;
template void save_wfn( wfn< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&, int) ;

template<typename s, int r, int c>
void load_wfn( wfn< s, r, c>& w, int cntl) {
/* 
  A controlling routine to read slater determinants from file.

  cntl : 
    0 - (default) read the first wfn. 
    N - Read N wfn Determinants 
*/
  std::ifstream F_IN ;

  if ( open_binary( F_IN) ) {
    qtzcntrl::shutdown( "IO Error (save_wfn)") ;
    }

  F_IN >> w.wfntyp ;
  F_IN >> w.e_scf ;
  read_eigen_bin( w.moc, F_IN) ;
  read_eigen_bin( w.eig, F_IN) ;

  F_IN.close() ;

  return ;

} ;

template void load_wfn( wfn< double, Eigen::Dynamic, Eigen::Dynamic>&, int) ;
template void load_wfn( wfn< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>&, int) ;

