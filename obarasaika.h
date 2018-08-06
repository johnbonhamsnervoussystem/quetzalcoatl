#include "basis.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#ifndef OBARASAIKA_H
#define OBARASAIKA_H

  double gauprm_ovl( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) ;

  double gauprm_T( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) ;

  double overlap_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) ;

  double kinetic_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) ;

  void ao_overlap( int n, basis_set& b) ;

  void ao_kinetic( int n, basis_set& b) ;

#endif
