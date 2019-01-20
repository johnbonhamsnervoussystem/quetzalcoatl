#include <Eigen/Core>
#include <queue>

#ifndef DIIS_H
#define DIIS_H

template<typename s>
class diis {
  private:
    /*
      n_e - Maximum Number of error vectors
      dim_e - dimension of error vectors
      w_diis - A queue of all the vectors
      B - A matrix for digonalization
      rho - The most recent density
      c - The diis vector for building an updated density
      v - The eigenvector of the diis
      t - temporary storage for incoming vectors
    */
    int n_e, dim_e ;
    std::vector<Eigen::Matrix< s, Eigen::Dynamic, 1>> w_diis ;
    Eigen::Matrix< s, Eigen::Dynamic, Eigen::Dynamic> B ;
    Eigen::Matrix< s, Eigen::Dynamic, 1> c, v, t, t1 ;

  public :
    diis( const int dim, const int n){
      dim_e = dim ;
      n_e = n ;
      w_diis.reserve( n_e) ;
      B.resize( n_e+1, n_e+1) ;
      c.resize( n_e+1) ;
      v.resize( n_e+1) ;
      t.resize( dim_e) ;
      t1.resize( dim_e) ;
      B.setZero() ;
      }

  void update( Eigen::Ref<Eigen::Matrix< s, Eigen::Dynamic, Eigen::Dynamic>> e, const int& extr) ;

} ;

#endif
