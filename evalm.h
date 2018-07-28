#include "constants.h"
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "common.h"
#include "hfwfn.h"
#include "tei.h"

#ifndef EVALM_H
#define EVALM_H

  int ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> G, const int nb) ;

  int ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> G, const int nb) ;

  int ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> pt, 
     Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> g, const int nb) ;

  int ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> pt, 
     Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> g, const int nb) ;

  int ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> G, const int nb) ;

  int ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> G, const int nb) ;

  int coulblk( std::vector<tei>& i, const Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> G, const int n) ;

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> G, const int nbasis) ;

  int coulblk( std::vector<tei>& i, const Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> G, const int n) ;

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> G, const int nbasis) ;

  double tranden ( common& c, hfwfn& a, hfwfn& b, Eigen::Ref<Eigen::MatrixXd> d) ;

  cd tranden ( common& c, hfwfn& a, hfwfn& b, Eigen::Ref<Eigen::MatrixXcd> d) ;

  double td_singleblock( int& o, int& n, Eigen::Ref<Eigen::MatrixXd> s1, Eigen::Ref<Eigen::MatrixXd> s2, Eigen::Ref<Eigen::MatrixXd> d) ;

  double obop ( common& c, Eigen::Ref<Eigen::MatrixXd> o, hfwfn& a, hfwfn& b) ;

  cd obop ( common& c, Eigen::Ref<Eigen::MatrixXcd> o, hfwfn& a, hfwfn& b) ;

  double fockop ( common& c, Eigen::Ref<Eigen::MatrixXd> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, double& O) ;

  cd fockop ( common& c, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, cd& O) ;

  cd fockop ( common& com, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, cd& ovl, std::ofstream& tfile, std::ofstream& ffile) ;

#endif
