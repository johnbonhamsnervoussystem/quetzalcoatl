#include "constants.h"
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "common.h"
#include "hfwfn.h"
#include "tei.h"

#ifndef EVALM_H
#define EVALM_H

  int ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> G, const int nb) ;

  int ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> G, const int nb) ;

  int ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> pt, 
     Eigen::Ref<Eigen::MatrixXf> p, Eigen::Ref<Eigen::MatrixXf> g, const int nb) ;

  int ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> pt, 
     Eigen::Ref<Eigen::MatrixXcf> p, Eigen::Ref<Eigen::MatrixXcf> g, const int nb) ;

  int ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> G, const int nb) ;

  int ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> G, const int nb) ;

  int coulblk( std::vector<tei>& i, const Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> G, const int n) ;

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> G, const int nbasis) ;

  int coulblk( std::vector<tei>& i, const Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> G, const int n) ;

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> G, const int nbasis) ;

  float tranden ( common& c, hfwfn& a, hfwfn& b, Eigen::Ref<Eigen::MatrixXf> d) ;

  cf tranden ( common& c, hfwfn& a, hfwfn& b, Eigen::Ref<Eigen::MatrixXcf> d) ;

  float obop ( common& c, Eigen::Ref<Eigen::MatrixXf> o, hfwfn& a, hfwfn& b) ;

  cf obop ( common& c, Eigen::Ref<Eigen::MatrixXcf> o, hfwfn& a, hfwfn& b) ;

  float fockop ( common& c, Eigen::Ref<Eigen::MatrixXf> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, float& O) ;

  cf fockop ( common& c, Eigen::Ref<Eigen::MatrixXcf> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, cf& O) ;

  cf fockop ( common& com, Eigen::Ref<Eigen::MatrixXcf> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, cf& ovl, std::ofstream& tfile, std::ofstream& ffile) ;

#endif
