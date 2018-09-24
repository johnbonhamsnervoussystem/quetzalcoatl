#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>

#ifndef SLADET_H
#define SLADET_H

template<typename s, int r, int c>
struct sladet {

  /* Track the type of wavefunction */
  int wfntyp ;
  double e_scf ;
  /* Store the wavefunction */
  Eigen::Matrix< s, r, c> moc ;
  Eigen::VectorXd eig ;

} ;

template<typename s, int r, int c>
void save_slater_det( sladet< s, r, c>& w, int cntl = 0) ;

template<typename s, int r, int c>
void write_to_bin( sladet< s, r, c>& w, std::ofstream& F_OUT) ;

template<typename s, int r, int c>
void load_slater_det( sladet< s, r, c>& w, int cntl = 0) ;

template<typename s, int r, int c>
void read_sladet( sladet< s, r, c>& w, std::ifstream& F_IN) ;

#endif
