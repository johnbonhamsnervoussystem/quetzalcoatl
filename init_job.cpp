#include "common.h"
#include <Eigen/Core>
#include "hubbard.h"
#include "init_job.h"
#include "nbodyint.h"
#include "qtzcntrl.h"
#include "pairing.h"
#include "r12.h"
#include "tei.h"
#include "util.h"
#include <vector>

template < class matrix>
void initialize( int wt, int ws, int hm, common& com, matrix& h, nbodyint<matrix>*& W, std::vector<tei>*& intarr, int& nbas) {
  matrix xs ;

  switch( hm){
    case 1 :
      /* 
        Molecular Hamiltonian 
      */
      if ( wt == 1){
        /* 
          Hartree-Fock
        */
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
          qtzcntrl::shutdown( " UHF NYI") ;
        } else {
          qtzcntrl::shutdown( " GHF NYI") ;
          }
      } else if (wt == 2){
        /* 
          Hartree-Fock-Bogoliubov
        */
        if ( ws == 1){
          h.resize( 2*nbas, 2*nbas) ;
          h.setZero() ;
          h.block( 0, 0, nbas, nbas).real() = com.getH() ;
          xs.resize( nbas, nbas) ;
          xs.real() = com.getXS() ;
          transform( 2, xs, h.block( 0, 0, nbas, nbas)) ;
          h.block( nbas, nbas, nbas, nbas) = -h.block( 0, 0, nbas, nbas) ;
          com.getr12( intarr) ;
          W = new r12<matrix>( intarr, xs, 4, nbas) ;
          xs.resize( 0, 0) ;
          break ;
        } else if ( ws == 2){
          qtzcntrl::shutdown( " UHFB NYI") ;
        } else {
          qtzcntrl::shutdown( " GHFB NYI") ;
          }
      } else {
        qtzcntrl::shutdown( " No other wavefunctions implemented ") ;
        }
    case 3 :
      /*
        Hubbard
      */
      if ( wt == 1){
        /* 
          Hartree-Fock
        */
        if ( ws == 1){
          h.resize( nbas, nbas) ;
          h.real() = com.getH() ;
          W = new hubbard<matrix>( com.getU(), 1, nbas) ;
          break ;
        } else {
          qtzcntrl::shutdown( " Only RHF Hubbard Implemented ") ;
          }
        }
        break ;
    case 4 :
      /* 
        Pairing Hamiltonian 
      */
      if ( wt == 1){
        /* 
          Hartree-Fock
        */
        if ( ws == 1){
          h.resize( nbas, nbas) ;
          h.real() = com.getH() ;
          W = new pairing<matrix>( com.getU(), 1, nbas) ;
        } else if ( ws == 2){
          qtzcntrl::shutdown( " UHF NYI") ;
        } else {
          qtzcntrl::shutdown( " GHF NYI") ;
          }
      } else if (wt == 2){
        /* 
          Hartree-Fock-Bogoliubov
        */
        if ( ws == 1){
          h.resize( 2*nbas, 2*nbas) ;
          h.setZero() ;
          h.block( 0, 0, nbas, nbas).real() = com.getH() ;
          h.block( nbas, nbas, nbas, nbas) = -h.block( 0, 0, nbas, nbas) ;
          W = new pairing<matrix>( com.getU(), 4, nbas) ;
        } else if ( ws == 2){
          qtzcntrl::shutdown( " UHF NYI") ;
        } else {
          h.resize( 2*nbas, 2*nbas) ;
          h.setZero() ;
          h.block( 0, 0, nbas, nbas).real() = com.getH() ;
          h.block( nbas, nbas, nbas, nbas) = h.block( 0, 0, nbas, nbas) ;
          W = new pairing<matrix>( com.getU(), 6, nbas) ;
          }
      } else {
        qtzcntrl::shutdown( " No other wavefunctions implemented ") ;
        }
      break ;
    default :
      qtzcntrl::shutdown( " Cannot initialize unknown Hamiltonian.") ;
      break ;
    }

  return ;

} ;

template void initialize( int, int, int, common&, Eigen::MatrixXd&, nbodyint<Eigen::MatrixXd>*&, std::vector<tei>*&, int&) ;

template void initialize( int, int, int, common&, Eigen::MatrixXcd&, nbodyint<Eigen::MatrixXcd>*&, std::vector<tei>*&, int&) ;
