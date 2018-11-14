#include "common.h"
#include <Eigen/Dense>
#include <iostream>
#include "wfn.h"
#include "tei.h"
#include "time_dbg.h"
#include <vector>

/*
  Correlated methods for post scf
*/

void mo_integrals( common &com, std::vector<tei>& intarr) {
/*
  Transform the integrals from the ao to the mo basis.
*/
  int nbasis = com.nbas() ;
  time_dbg mo_integrals_time = time_dbg("mo_integrals") ;
/*
  Start with real integrals
*/
  wfn< double, Eigen::Dynamic, Eigen::Dynamic> w ;

  w.moc.resize( nbasis, nbasis) ;
  w.eig.resize( nbasis) ;

  load_wfn( w) ;

  std::cout << w.moc << std::endl ;

  mo_integrals_time.end() ;

  return ;

} /* End mo_integrals */

