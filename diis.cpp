#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include "diis.h"

template<typename s>
void diis<s>::update( Eigen::Ref<Eigen::Matrix< s, Eigen::Dynamic, Eigen::Dynamic>> p, const int& extr){

  /*
    If extr == 1 then update the density using DIIS.  Otherwise we just add it to the
      queue.
  */

  int i, j, np ;
  /*
    Convert the incoming density to linear. It's column major by default.
  */
  t = Eigen::Map<Eigen::Matrix< s, Eigen::Dynamic, 1>>( p.data(), p.cols()*p.rows(), 1) ;
  
  /*
    Add the density vector to the queue
  */
  w_diis.push_back( t) ;
  np = w_diis.size() ;

  /*
    Do the DIIS 
  */
  if ( np > 2 && extr == 1){
    B.block(0, 0, np, np).setConstant( 1.0e0) ;
    B( np-1, np-1) = 0.0e0 ;
    for ( i = 0; i < np - 1; i++){
      t = w_diis[i+1].conjugate() - w_diis[i].conjugate() ;
      for ( j = i; j < np - 1; j++){
        t1 = w_diis[j+1] - w_diis[j] ;
        B( i, j) = t.dot( t1) ;
        B( j, i) = std::conj(B( i, j)) ;
        }
      }
    v.setZero() ;
    c.setZero() ;
    v( np-1) = 1.0e0 ;
    c.head( np) = B.block( 0, 0, np, np).colPivHouseholderQr().solve(v.head( np)) ;

    t.setZero() ;

    for( i = 0; i < np-1; i++){
      t += c(i)*w_diis[i] ;
      }

    /*
      t is our new density send it back into p
    */
    p = Eigen::Map<Eigen::Matrix< s, Eigen::Dynamic, Eigen::Dynamic>>( t.data(), p.rows(), p.cols()) ;
    }

  while ( np > n_e ){
    w_diis.erase( w_diis.begin()) ;
    np = w_diis.size() ;
    }

  return ;

} ;

template void diis<double>::update( Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>, const int&) ;
template void diis<std::complex<double>>::update( Eigen::Ref<Eigen::Matrix< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>, const int&) ;

