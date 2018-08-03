#include "basis.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#ifndef OBARASAIKA_H
#define OBARASAIKA_H

  double gauprm_ovl( const double& za, const double& ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, const double& zb, const double& cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) ;

  double overlap_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) ;

  void ao_overlap( int n, basis_set& b) ;

#endif
