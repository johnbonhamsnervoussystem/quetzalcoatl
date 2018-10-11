#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>

#ifndef WFN_H
#define WFN_H

//template<typename s, int r, int c>
//struct sladet {
//
//  /* Track the type of wavefunction */
//  int wfntyp ;
//  double e_scf ;
//  /* Store the wavefunction */
//  Eigen::Matrix< s, r, c> moc ;
//  Eigen::VectorXd eig ;
//
//} ; 

template<typename s, int r, int c>
struct wfn {

  /* Track the type of wavefunction */
  int wfntyp ;
  double e_scf ;
  /* Store the wavefunction */
  Eigen::Matrix< s, r, c> moc ;
  Eigen::VectorXd eig ;

} ;

template<typename s, int r, int c>
//void save_slater_det( sladet< s, r, c>& w, int cntl = 0) ;
void save_slater_det( wfn< s, r, c>& w, int cntl = 0) ;

template<typename s, int r, int c>
//void load_slater_det( sladet< s, r, c>& w, int cntl = 0) ;
void load_slater_det( wfn< s, r, c>& w, int cntl = 0) ;

#endif
