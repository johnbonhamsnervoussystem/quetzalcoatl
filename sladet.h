#include <Eigen/Core>
#include <Eigen/Dense>

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
void write_to_bin( sladet<s, r, c>& w) ;

#endif
