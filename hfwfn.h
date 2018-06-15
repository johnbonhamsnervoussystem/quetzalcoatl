#include <Eigen/Dense>
#include <vector>
#include <string>
#include "common.h"
#include "tei.h"

#ifndef HFWFN_H
#define HFWFN_H

class hfwfn {
  private :
    float energy ;
    int wfn_ityp ;
    std::string wfn_styp ;
    Eigen::VectorXf eig_v ;
    Eigen::MatrixXf mo_rcof ;
    Eigen::MatrixXcf mo_ccof ;

  public :

  void init ( common& com, std::vector<tei>& intarr, std::string wfn) ;
  void fil_mos ( int n, Eigen::Ref<Eigen::MatrixXf> mo, int wfn) ;
  void fil_mos ( int n, Eigen::Ref<Eigen::MatrixXcf> mo, int wfn) ;

  void set_mos ( Eigen::Ref<Eigen::MatrixXf> mo) ;
  void set_mos ( Eigen::Ref<Eigen::MatrixXcf> mo) ;

  void get_mos ( Eigen::Ref<Eigen::MatrixXf> mo) ;
  void get_mos ( Eigen::Ref<Eigen::MatrixXcf> mo) ;

/*
 *  1 -> rrhf
 *  2 -> crhf
 *  3 -> ruhf
 *  4 -> cuhf
 *  5 -> rghf
 *  6 -> cghf
 *           */

  int get_wti ( void) ;
  std::string get_wts ( void) ;

  void prt_mos( void) ;
  void prt_eig( void) ;
  void prt_ene( float nn = 0.0) ;

void ia ( Eigen::Ref<Eigen::MatrixXf> mo, int i, int a) ;

void ia ( Eigen::Ref<Eigen::MatrixXcf> mo, int i, int a) ;

} ;

#endif
