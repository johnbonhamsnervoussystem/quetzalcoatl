#include <Eigen/Core>
#include <queue>

#ifndef DIIS_H
#define DIIS_H

template<typename s>
class diis {
/*
  Based on the papers :
    Pulay, Péter. "Convergence acceleration of iterative sequences. 
    The case of SCF iteration." Chemical Physics Letters 73, no. 2 
    (1980): 393-398.

    Pulay, Péter. "Improved SCF convergence aceleration." 
    Journal of Computational Chemistry 3(4), (1982): 556-560.
*/
  private:
    /*
      n_e - Maximum Number of error vectors
      dim_e - dimension of error vectors
      e_diis - A queue of all the error vectors
      e_diis - A queue of all the density vectors
      B - A matrix for digonalization
      rho - The most recent density
      c - The diis vector for building an updated density
      v - The eigenvector of the diis
      t - temporary storage for incoming vectors
    */
    int n_e, dim_e ;
    bool diis_on ;
    std::vector<Eigen::Matrix< s, Eigen::Dynamic, 1>> e_diis, f_diis ;
    Eigen::Matrix< s, Eigen::Dynamic, Eigen::Dynamic> B, scr ;
    Eigen::Matrix< s, Eigen::Dynamic, 1> c, v, t, t1 ;

  public :
    diis( const int dim, const int n, const int type){
      diis_on = false ;
      if ( type == 1 ){
      /* Real diis */
        dim_e = dim*dim ;
        n_e = n ;
        e_diis.reserve( n_e) ;
        f_diis.reserve( n_e) ;
        B.resize( n_e+1, n_e+1) ;
        c.resize( n_e+1) ;
        v.resize( n_e+1) ;
        t.resize( dim_e) ;
        t1.resize( dim_e) ;
        scr.resize( dim, dim) ;
      } else {
      /* complex diis */
        dim_e = dim*dim ;
        n_e = n ;
        e_diis.reserve( n_e) ;
        f_diis.reserve( n_e) ;
        B.resize( 2*n_e+2, 2*n_e+2) ;
        c.resize( 2*n_e+2) ;
        v.resize( 2*n_e+2) ;
        t.resize( dim_e) ;
        t1.resize( dim_e) ;
        scr.resize( dim, dim) ;
        }
      }

  void update( Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> f, int& d_cntrl) ;
  void update( Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> f, int& d_cntrl) ;

} ;

class diis_control{
/*
  A container to hold various control options for passing into the main routines and declutter the 
  function calls.
*/
  private :
    /* What are we using to check the diis*/
    int control_switch = 0 ;

  public :
    /* Are we doing diis in the routine? */
    bool do_diis ;
    /* Is diis currently active? */
    int diis_switch ;
    /* Number of diis vectors to use */
    int ndiisv ;
    /* Real or complex diis */
    int diistype ;

    /* Threshold value for the energy difference test */
    double edif_v ;

    /* Value to check the energy difference against */
    double cntl_dif ;
    /* Value to chekc the threshold against */
    double cntl_ref ;

    void set_switch( int opt) ;
    void toggle( double e) ;

} ;

#endif
