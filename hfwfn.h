#include <Eigen/Dense>
#include <vector>
#include <string>
#include "common.h"
#include "tei.h"

#ifndef HFWFN_H
#define HFWFN_H

class hfwfn {
  private :
    double energy ;
    int wfn_ityp ;
    std::string wfn_styp ;
    Eigen::VectorXd eig_v ;
    Eigen::MatrixXd mo_rcof ;
    Eigen::MatrixXcd mo_ccof ;

  public :

//  void init ( common& com, std::vector<tei>& intarr, std::string wfn) ;
  void fil_mos ( int n, Eigen::Ref<Eigen::MatrixXd> mo, int wfn) ;
  void fil_mos ( int n, Eigen::Ref<Eigen::MatrixXcd> mo, int wfn) ;

  void set_mos ( Eigen::Ref<Eigen::MatrixXd> mo) ;
  void set_mos ( Eigen::Ref<Eigen::MatrixXcd> mo) ;

  void get_mos ( Eigen::Ref<Eigen::MatrixXd> mo) ;
  void get_mos ( Eigen::Ref<Eigen::MatrixXcd> mo) ;

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
  void prt_ene( double nn = 0.0) ;

void ia ( Eigen::Ref<Eigen::MatrixXd> mo, int i, int a) ;

void ia ( Eigen::Ref<Eigen::MatrixXcd> mo, int i, int a) ;

} ;

#endif
