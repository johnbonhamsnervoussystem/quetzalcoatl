#include "common.h"
#include <Eigen/Core>
#include "init_job.h"
#include "nbodyint.h"
#include "qtzcntrl.h"
#include "pairing.h"
#include "r12.h"
#include "tei.h"
#include "util.h"
#include <vector>

template < class matrix>
void initialize( int ws, int hm, common& com, matrix& h, nbodyint<matrix>*& W, std::vector<tei>*& intarr, int& nbas) {
  matrix xs ;

  switch( hm){
    case 1 :
      /* Molecular Hamiltonian */
      if ( ws == 1){
        h.resize( nbas, nbas) ;
        h.real() = com.getH() ;
        xs.resize( nbas, nbas) ;
        xs.real() = com.getXS() ;
        transform( 2, xs, h) ;
        com.getr12( intarr) ;
        W = new r12<matrix>( intarr, xs, 1, nbas) ;
        xs.resize( 0, 0) ;
        break ;
      } else if ( ws == 2){
      } else {
        }
    case 4 :
        h.resize( nbas, nbas) ;
        h.real() = com.getH() ;
        W = new pairing<matrix>( com.getU(), 1, nbas) ;
      break ;
    default :
      qtzcntrl::shutdown( " Cannot initialize unknown Hamiltonian.") ;
      break ;
    }

  return ;

} ;

template void initialize( int, int, common&, Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*&, std::vector<tei>*&, int&) ;

template void initialize( int, int, common&, Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*&, std::vector<tei>*&, int&) ;
