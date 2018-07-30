#include "basis.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "obarasaika.h"
#include <string>
#include <vector>

  basis_set load_sto1g ( Eigen::VectorXi AtN, Eigen::MatrixXd c) {
    /* STO-3G basis set
     cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
           doi: 10.1063/1.1672392
     obtained from https://bse.pnl.gov/bse/portal */
    gau_prm g ;
    sto s ;
    double n ;
    double t3_f4 = 3.0e0/4.0e0 ;
    basis_fc bf ;
    basis_set basis ;
    basis.nbas = 0 ;

    for( int a=0; a<AtN.size(); ++a) {

      switch (AtN(a)) {
        case 1: // Z=1: hydrogen
          /* One s orbital */
          g.x = 0.270950 ;
          n = (2*g.x)/pi ;
          g.c = d1*pow( n, t3_f4) ;
          s.g.push_back(g) ;
          s.g.push_back(g) ;
          s.l.setZero() ;

          s.nprm = 1 ;
          bf.s.push_back(s) ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
 
          basis.b.push_back(bf) ;
          basis.nbas ++ ;

        }
      }

  return basis ;

 } ;

  basis_set load_sto3g ( Eigen::VectorXi AtN, Eigen::MatrixXd c) {
    /* STO-3G basis set
     cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
           doi: 10.1063/1.1672392
     obtained from https://bse.pnl.gov/bse/portal */
    gau_prm g ;
    sto s ;
    double n ;
    double t3_f4 = 3.0e0/4.0e0 ;
    basis_fc bf ;
    basis_set basis ;
    basis.nbas = 0 ;


    for( int a=0; a<AtN.size(); ++a) {

      switch (AtN(a)) {
        case 1: // Z=1: hydrogen
          /* One s orbital */
          g.x = 3.425250914 ;
          n = (2*g.x)/pi ;
          g.c = 0.1543289673*pow( n, t3_f4) ;
          s.g.push_back(g) ;
          g.x = 0.6239137298 ;
          n = (2*g.x)/pi ;
          g.c = 0.5353281423*pow( n, t3_f4) ;
          s.g.push_back(g) ;
          g.x = 0.1688554040 ;
          n = (2*g.x)/pi ;
          g.c = 0.4446345422*pow( n, t3_f4) ;
          s.g.push_back(g) ;
          s.l.setZero() ;

          s.nprm = 3 ;
          bf.s.push_back(s) ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
 
          basis.b.push_back(bf) ;
          basis.nbas ++ ;

        }
      }

  return basis ;

 } ;

  basis_set build_basis ( std::string bas_name, Eigen::Ref<Eigen::VectorXi> AtN, Eigen::Ref<Eigen::MatrixXd> coord) {
  /* Load a basis set */
  const std::string b_sto3g = "sto-3g" ;
  const std::string b_sto1g = "sto-1g" ;
  basis_set basis ;

  if ( bas_name == b_sto3g ){
    basis = load_sto3g( AtN, coord) ;
    }

  if ( bas_name == b_sto1g ){
    basis = load_sto1g( AtN, coord) ;
    }

   norm_basis ( basis) ;

   return basis ;

   } ;

  void norm_basis ( basis_set& b) {
  /* Normalize the basis */
    int nsto ;
    int nbas = b.nbas ;
    double s = 0.0 ;

    for( int i = 0; i < nbas; i++ ){
      nsto = b.b[i].nshl ;
      for( int j = 0; j < nsto; j++){
        s = overlap_sto( b.b[i].s[j], b.b[i].c, b.b[i].s[j], b.b[i].c) ;
        b.b[i].s[j].norm = d1/sqrt(s) ;
	}
      } ;

    } ;

  void print_basis( basis_set& basis){
    /* Dump basis information for debugging */
    int nprm ;
    int nsto ;
    int nbas = basis.nbas ;

    std::cout << "basis = " << nbas << std::endl ;
    for( int i = 0; i < nbas; i++ ){
      std::cout << "basis = " << i << std::endl ;
      nsto = basis.b[i].nshl ;
      for( int j = 0; j < nsto; j++){
        nprm = basis.b[i].s[j].nprm ;
        std::cout << "lx = " << basis.b[i].s[j].l(0) << std::endl ;
        std::cout << "ly = " << basis.b[i].s[j].l(1) << std::endl ;
        std::cout << "lz = " << basis.b[i].s[j].l(2) << std::endl ;
        std::cout << " cof    " << " exp    " << std::endl ;
        std::cout << "----------------" << std::endl ;
        for( int k = 0; k < nprm; k++){
          std::cout << basis.b[i].s[j].g[k].c << " " ;
          std::cout << basis.b[i].s[j].g[k].x << std::endl ;
          }
        }

      } ;

    return ;

    }

