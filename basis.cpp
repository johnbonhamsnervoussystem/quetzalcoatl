#include "basis.h"
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>

  basis_set load_sto3g ( Eigen::VectorXi AtN, Eigen::MatrixXd c) {
    /* STO-3G basis set
     cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
           doi: 10.1063/1.1672392
     obtained from https://bse.pnl.gov/bse/portal */
    gau_prm g ;
    sto s ;
    basis_fc bf ;
    basis_set basis ;
    basis.nbas = 0 ;

    for( int a=0; a<AtN.size(); ++a) {

      switch (AtN(a)) {
        case 1: // Z=1: hydrogen
          /* One s orbital */
          g.x = 3.425250910 ;
          g.c = 0.15432897 ;
          s.g.push_back(g) ;
          g.x = 0.623913730 ;
          g.c = 0.53532814 ;
          s.g.push_back(g) ;
          g.x = 0.168855400 ;
          g.c = 0.44463454 ;
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
  basis_set basis ;

  if ( bas_name == b_sto3g ){
    basis = load_sto3g( AtN, coord) ;
    }

   return basis ;

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

