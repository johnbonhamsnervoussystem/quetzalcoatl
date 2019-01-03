#include "basis.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "tei.h"
#include <vector>

#ifndef OBARASAIKA_H
#define OBARASAIKA_H

  double gauprm_ovl( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) ;

  double gauprm_T( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) ;

  double gauprm_V( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb, Eigen::Ref<Eigen::Vector3d> c, double q ) ; 

  double gauprm_r12( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb, double zc, double cc, Eigen::Ref<Eigen::Vector3d> c, Eigen::Ref<Eigen::Vector3i> lc, double zd, double cd, Eigen::Ref<Eigen::Vector3d> d, Eigen::Ref<Eigen::Vector3i> ld )  ;

  double overlap_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) ;

  double kinetic_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) ;

  double nucelecV_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb, Eigen::Ref<Eigen::MatrixXd> n_c, Eigen::Ref<Eigen::VectorXd> q) ;

  void ao_overlap( int n, basis_set& b, Eigen::Ref<Eigen::MatrixXd> ovl) ;

  void ao_kinetic( int n, basis_set& b, Eigen::Ref<Eigen::MatrixXd> T) ;

  void ao_eN_V( int natm, basis_set& b, Eigen::Ref<Eigen::MatrixXd> n_c, Eigen::Ref<Eigen::VectorXd> q, Eigen::Ref<Eigen::MatrixXd> V) ; 

/*  Eigen::Tensor< double, 4> ao_tei( int natm, basis_set& b ) ;
  void ao_tei( int nbas, int natm, basis_set& b, Eigen::Ref<Eigen::Tensor< double, nbas>> eri ) ; */

  void list_ao_tei( int natm, basis_set& b, std::vector<tei>& intarr) ;

#endif
