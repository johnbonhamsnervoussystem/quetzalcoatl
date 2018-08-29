/* Primary routine for Quetzalcoatl */
#include "basis.h"
#include "constants.h"
#include <complex>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <sstream>
#include <fstream>
#include "binio.h"
#include "common.h"
#include "evalm.h"
#include "hfrout.h"
#include "integr.h"
#include "obarasaika.h"
#include "project.h"
#include "qtzio.h"
#include "solver.h"
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include "wigner.h"

int main(int argc, char *argv[]) {
/*
 * This is my own implementation of various electronic stucture methods.  
 * There is lots of work to do and improvements to be made and the purpose
 * is supposed to be minamly pedagoical as well as producing results. 
 * Currently everything is done in memory.
 * To Do :: 
 *   - Parse options and set route through program.
 *   -implement stability 
 *   -matrix elements between determinants
 *   -Eventually, I will need to write a routine which grabs a chunk of memory,
 *   rather than constantly allocating and dealllocating.  For now, I need things
 *   that work and produce useful results.  
 *   - Debug spin rotations on a determinant.
 *   - Interface with non-gaussian routines for matrix elements
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
  common com ;
  std::vector<tei> intarr ;

  int job=0  ;
  int i, j, k, l, nbas ;
  double cx = 0.0e0 ;
  double cy = 0.0e0 ;
  double cz = 0.0e0 ;
  double r = 0.0e0 ;
  double r2 = 0.0e0 ;
  double n_rep = 0.0e0 ;
  double val = 0.0e0 ;
  cd ejunk ;
  cd ojunk ;
  hfwfn det1 ;
  hfwfn det2 ;
  std::vector<hfwfn> tst_vec ;
  std::vector<std::string> wfn_vec ;
  std::string trden="tden.rwf" ;
  std::string fokmat="fmat.rwf" ;
  Eigen::MatrixXd c ;
  Eigen::VectorXd a ;
  Eigen::MatrixXd S ;
  Eigen::MatrixXd T ;
  Eigen::MatrixXd V ;
  Eigen::Tensor< double, 4> eri( 0, 0, 0, 0) ;
  basis_set b ; 
  time_dbg quetz_time = time_dbg("Quetzalcoatl") ;

  /* File reading and header variables. */
  std::stringstream ss ;
  std::string inpfile ;

  ss << argv[1] ;
  ss >> inpfile ;
  read_input( com, inpfile) ;

  std::cout << " ---   ---   ---   --- " << std::endl << "|" << std::endl 
    << "|  Quetzalcoatl v1.0  " << std::endl ;
  std::cout << "| " << std::endl ;
  std::cout << " ---   ---   ---   --- " << std::endl << "|" << std::endl ;
  std::cout << "|  Reading from input : " << inpfile << std::endl ;
  std::cout << "|  Basis Set : " << com.bnam() << std::endl ;
  std::cout << "| " << std::endl ;

  std::cout << "| " << std::endl ;
  std::cout << " ---   ---   ---   --- " << std::endl ;

  /* Step 1 :
     Build the relevant data in memory.
     SCF routines
     real/complex reastricted
                  unrestricted
                  generalized  */
  

  c.resize( com.natm(), 3) ;
  a.resize( com.natm()) ;
  c = com.getC() ;
  a = com.getA() ;
 
//  com.nrep( n_rep ) ;

  b = build_basis( com.bnam(), a, c) ;
  nbas = b.nbas ;
  com.nbas( nbas) ;
  std::cout << "nbas " << nbas << std::endl ;
  S.resize( nbas, nbas) ;
  ao_overlap( com.natm(), b, S) ;
  T.resize( nbas, nbas) ;
  ao_kinetic( com.natm(), b, T) ;
  V.resize( nbas, nbas) ;
  ao_eN_V( com.natm(), b, c, a, V) ;
  list_ao_tei( com.natm(), b, intarr) ; 

  int opt = 1 ;
  scf_drv( com, T, V, S, intarr, opt) ;

  /* Deallocate the memory and exit. */
  intarr.clear() ;
  V.resize( 0, 0) ; 
  T.resize( 0, 0) ; 
  S.resize( 0, 0) ; 
  quetz_time.end() ;

  return 0 ;

} /* End Quetzacoatl */

