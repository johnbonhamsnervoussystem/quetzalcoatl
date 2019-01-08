#include "constants.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "common.h"
#include "tei.h"

#ifndef EVALM_H
#define EVALM_H

void ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> G, const int nb) ;

void ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> p, 
  Eigen::Ref<Eigen::MatrixXcd> G, const int nb) ;

void ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> pt, 
   Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> g, const int nb) ;

void ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> pt, 
   Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> g, const int nb) ;

void ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> G, const int nb) ;

void ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> p, 
  Eigen::Ref<Eigen::MatrixXcd> G, const int nb) ;

void ctrPairg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> k, 
  Eigen::Ref<Eigen::MatrixXd> D, const int nb) ;

void ctrPairg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> k, 
  Eigen::Ref<Eigen::MatrixXcd> D, const int nb) ;

int coulblk( std::vector<tei>& i, const Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> G, const int n) ;

int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> G, const int nbasis) ;

int coulblk( std::vector<tei>& i, const Eigen::Ref<Eigen::MatrixXcd> p, 
  Eigen::Ref<Eigen::MatrixXcd> G, const int n) ;

int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
  Eigen::Ref<Eigen::MatrixXcd> G, const int nbasis) ;

int pairing_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> D, const int nbasis) ;

int pairing_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
  Eigen::Ref<Eigen::MatrixXcd> D, const int nbasis) ;

int twin_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> D, const int nbasis) ;

int twin_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
  Eigen::Ref<Eigen::MatrixXcd> D, const int nbasis) ;

  double tranden1 ( int& ne, int& nb, Eigen::Ref<Eigen::MatrixXd> wa, Eigen::Ref<Eigen::MatrixXd> wb, Eigen::Ref<Eigen::MatrixXd> dmat, Eigen::Ref<Eigen::MatrixXd> s) ;

  cd tranden1 ( int& ne, int& nb, Eigen::Ref<Eigen::MatrixXcd> wa, Eigen::Ref<Eigen::MatrixXcd> wb, Eigen::Ref<Eigen::MatrixXcd> dmat, Eigen::Ref<Eigen::MatrixXcd> s) ;

double tranden2 ( int& ne1, int& ne2, int& nbasis, Eigen::Ref<Eigen::MatrixXd> wfna, Eigen::Ref<Eigen::MatrixXd> wfnb, Eigen::Ref<Eigen::MatrixXd> dabmat, Eigen::Ref<Eigen::MatrixXd> scr1) ;

cd tranden2 ( int& nele, int& nbasis, Eigen::Ref<Eigen::MatrixXcd> wfna, Eigen::Ref<Eigen::MatrixXcd> wfnb, Eigen::Ref<Eigen::MatrixXcd> dabmat, Eigen::Ref<Eigen::MatrixXcd> scr) ;

//double obop ( common& c, Eigen::Ref<Eigen::MatrixXd> o, hfwfn& a, hfwfn& b) ;

//cd obop ( common& c, Eigen::Ref<Eigen::MatrixXcd> o, hfwfn& a, hfwfn& b) ;

//double fockop ( common& c, Eigen::Ref<Eigen::MatrixXd> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, double& O) ;

//cd fockop ( common& c, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, cd& O) ;

//cd fockop ( common& com, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, hfwfn& a, hfwfn& b, cd& ovl, std::ofstream& tfile, std::ofstream& ffile) ;

#endif
