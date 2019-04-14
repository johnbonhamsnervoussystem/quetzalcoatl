#include "basis.h"
#include "common.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "obarasaika.h"
#include <string>
#include "util.h"
#include <vector>

  basis_set load_sto1g ( Eigen::VectorXd AtN, Eigen::MatrixXd c) {
    /* STO-3G basis set
     cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
           doi: 10.1063/1.1672392
     obtained from https://bse.pnl.gov/bse/portal */
    gau_prm g ;
    sto s ;
    atm_basis bf ;
    basis_set basis ;
    basis.nbas = 0 ;

    for( int a=0; a<AtN.size(); ++a) {

      bf.s.clear() ;
      switch ( static_cast<int>(AtN(a)+0.5)) {
        case 1: // Z=1: hydrogen
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 0.270950 ;
          g.c = norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;

          bf.s.push_back(s) ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
 
          basis.b.push_back(bf) ;

          basis.nbas++ ;

          break ;

        case 6: // Z=6: Carbon
          /* 2px orbital */
          s.g.clear() ;
          s.l << 1,0,0 ;
          g.x = 2.9412494 ;
          g.c = norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
          basis.b.push_back(bf) ;
          break ;

        case 8: // Z=8: Oxygen
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 130.7093200 ;
          g.c = 0.15432897*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 23.8088610 ;
          g.c = 0.53532814*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 6.4436083 ;
          g.c = 0.44463454*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
          basis.b.push_back(bf) ;
          break ;

        default :
          std::cout << "Atomic number " << AtN(a) ; 
          std::cout << " not currently supported."  << std::endl ;
          break ;

        }
      }

  return basis ;

 } ;

  basis_set load_sto3g ( Eigen::VectorXd AtN, Eigen::MatrixXd c) {
    /* STO-3G basis set
     cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
           doi: 10.1063/1.1672392
     obtained from https://bse.pnl.gov/bse/portal */
    gau_prm g ;
    sto s ;
    atm_basis bf ;
    basis_set basis ;
    basis.nbas = 0 ;


    for( int a=0; a<AtN.size(); a++) {

      switch ( static_cast<int>(AtN(a)+0.5)) {
        case 1: // Z=1: hydrogen
          /* One s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 3.425250914 ;
          g.c = 0.1543289673*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.6239137298 ;
          g.c = 0.5353281423*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.1688554040 ;
          g.c = 0.4446345422*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
          basis.b.push_back(bf) ;
          basis.nbas++ ;
          break ;

        case 2: // Z=2: Helium
          /* One s orbital */
          s.l.setZero() ;
          s.g.clear() ;
          bf.s.clear() ;
          g.x = 6.36242139 ;
          g.c = 0.15432897*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.15892300 ;
          g.c = 0.53532814*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.31364979 ;
          g.c = 0.44463454*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;

          bf.s.push_back(s) ;

          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 1 ;
 
          basis.b.push_back(bf) ;
          break ;

        case 6: // Z=6: Carbon
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 71.6168370 ;
          g.c = 0.15432897*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 13.0450960 ;
          g.c = 0.53532814*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 3.5305122 ;
          g.c = 0.44463454*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          /* 2s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 2.9412494 ;
          g.c = -0.09996723*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.6834831 ;
          g.c = 0.39951283*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2222899 ;
          g.c = 0.70011547*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2px orbital */
          s.g.clear() ;
          s.l << 1,0,0 ;
          g.x = 2.9412494 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.6834831 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2222899 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2py orbital */
          s.g.clear() ;
          s.l << 0,1,0 ;
          g.x = 2.9412494 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.6834831 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2222899 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2pz orbital */
          s.g.clear() ;
          s.l << 0,0,1 ;
          g.x = 2.9412494 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.6834831 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2222899 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 5 ;
          basis.b.push_back(bf) ;
          break ;

        case 7: // Z=7: Nitrogen
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 99.1061690 ;
          g.c = 0.15432897*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 18.0523120 ;
          g.c = 0.53532814*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 4.8856602 ;
          g.c = 0.44463454*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          /* 2s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 3.7804559 ;
          g.c = -0.09996723*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.8784966 ;
          g.c = 0.39951283*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2857144 ;
          g.c = 0.70011547*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2px orbital */
          s.g.clear() ;
          s.l << 1,0,0 ;
          g.x = 3.7804559 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.8784966 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2857144 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2py orbital */
          s.g.clear() ;
          s.l << 0,1,0 ;
          g.x = 3.7804559 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.8784966 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2857144 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2pz orbital */
          s.g.clear() ;
          s.l << 0,0,1 ;
          g.x = 3.7804559 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.8784966 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.2857144 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 5 ;
          basis.b.push_back(bf) ;
          break ;

        case 8: // Z=8: Oxygen
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 130.7093200 ;
          g.c = 0.15432897*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 23.8088610 ;
          g.c = 0.53532814*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 6.4436083 ;
          g.c = 0.44463454*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          /* 2s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 5.0331513 ;
          g.c = -0.09996723*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.1695961 ;
          g.c = 0.39951283*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.3803890 ;
          g.c = 0.70011547*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2px orbital */
          s.g.clear() ;
          s.l << 1,0,0 ;
          g.x = 5.0331513 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.1695961 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.3803890 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2py orbital */
          s.g.clear() ;
          s.l << 0,1,0 ;
          g.x = 5.0331513 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.1695961 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.3803890 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2pz orbital */
          s.g.clear() ;
          s.l << 0,0,1 ;
          g.x = 5.0331513 ;
          g.c = 0.15591627*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.1695961 ;
          g.c = 0.60768372*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.3803890 ;
          g.c = 0.39195739*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 5 ;
          basis.b.push_back(bf) ;
          break ;

        default :
          std::cout << "Atomic number " << AtN(a) ; 
          std::cout << " not currently supported."  << std::endl ;
          break ;

        }
      }

  return basis ;

 } ;

  basis_set load_cc_pVDZ ( Eigen::VectorXd AtN, Eigen::MatrixXd c) {
/*
  cc-pVDZ  EMSL  Basis Set Exchange Library   12/16/18 2:18 PM
  Elements                             References
  --------                             ----------
  H      : T.H. Dunning, Jr. J. Chem. Phys. 90, 1007 (1989).
  He     : D.E. Woon and T.H. Dunning, Jr. J. Chem. Phys. 100, 2975 (1994).
  Li - Ne: T.H. Dunning, Jr. J. Chem. Phys. 90, 1007 (1989).
  Na - Mg: D.E. Woon and T.H. Dunning, Jr.  (to be published)
  Al - Ar: D.E. Woon and T.H. Dunning, Jr.  J. Chem. Phys. 98, 1358 (1993).
  Sc - Zn: N.B. Balabanov and K.A. Peterson, J. Chem. Phys. 123, 064107 (2005),
  N.B. Balabanov and K.A. Peterson, J. Chem. Phys. 125, 074110 (2006)
  Ca     : J. Koput and K.A. Peterson, J. Phys. Chem. A, 106, 9595 (2002).
*/
    gau_prm g ;
    sto s ;
    atm_basis bf ;
    basis_set basis ;
    basis.nbas = 0 ;


    for( int a=0; a<AtN.size(); a++) {

      switch ( static_cast<int>(AtN(a)+0.5)) {

        case 2: // Z=2: Helium
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 38.36e0 ;
          g.c = 0.0238090*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 5.77e0 ;
          g.c = 0.1548910*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.24 ;
          g.c = 0.4699870*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 0.2976 ;
          g.c = d1*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2px orbital */
          s.g.clear() ;
          s.l << 1,0,0 ;
          g.x = 1.275 ;
          g.c = d1*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2py orbital */
          s.g.clear() ;
          s.l << 0,1,0 ;
          g.x = 1.275 ;
          g.c = d1*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          /* 2pz orbital */
          s.g.clear() ;
          s.l << 0,0,1 ;
          g.x = 1.275 ;
          g.c = d1*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 5 ;
          basis.b.push_back(bf) ;
          break ;

        default :
          std::cout << "Atomic number " << AtN(a) ; 
          std::cout << " not currently supported."<< std::endl ;
          break ;

        }
      }

  return basis ;

 } ;

  basis_set load_6_31g ( Eigen::VectorXd AtN, Eigen::MatrixXd c) {
/*
  6-31G  EMSL  Basis Set Exchange Library   12/16/18 3:19 PM
  Elements                             References
  --------                             ----------
  H - He: W.J. Hehre, R. Ditchfield and J.A. Pople, J. Chem. Phys. 56,
  Li - Ne: 2257 (1972).  Note: Li and B come from J.D. Dill and J.A.
  Pople, J. Chem. Phys. 62, 2921 (1975).
  Na - Ar: M.M. Francl, W.J. Petro, W.J. Hehre, J.S. Binkley, M.S. Gordon,
  D.J. DeFrees and J.A. Pople, J. Chem. Phys. 77, 3654 (1982)
  K  - Zn: V. Rassolov, J.A. Pople, M. Ratner and T.L. Windus, J. Chem. Phys.
  109, 1223 (1998)
  Note: He and Ne are unpublished basis sets taken from the Gaussian
  program
*/
    gau_prm g ;
    sto s ;
    atm_basis bf ;
    basis_set basis ;
    basis.nbas = 0 ;


    for( int a=0; a<AtN.size(); a++) {

      switch ( static_cast<int>(AtN(a)+0.5)) {

        case 1: // Z=1: Hydrogen
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 18.7311370e0 ;
          g.c = 0.03349460*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 2.8253937e0 ;
          g.c = 0.23472695*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 0.6401217e0 ;
          g.c = 0.81375733*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 0.1612778e0 ;
          g.c = d1*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 2 ;
          basis.b.push_back(bf) ;
          break ;

        case 2: // Z=2: Helium
          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          bf.s.clear() ;
          g.x = 38.4216340 ;
          g.c = 0.04013973935*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 5.7780300 ;
          g.c = 0.2612460970*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          g.x = 1.2417740 ;
          g.c = 0.7931846246*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 3 ;
          bf.s.push_back(s) ;

          basis.nbas++ ;

          /* 1s orbital */
          s.g.clear() ;
          s.l.setZero() ;
          g.x = 0.2979640 ;
          g.c = d1*norm_gprm( g.x, s.l) ;
          s.g.push_back(g) ;
          s.nprm = 1 ;
          bf.s.push_back(s) ;
 
          basis.nbas++ ;

          bf.c( 0) = c( a, 0) ;
          bf.c( 1) = c( a, 1) ;
          bf.c( 2) = c( a, 2) ;
          bf.nshl = 2 ;
          basis.b.push_back(bf) ;
          break ;

        default :
          std::cout << "Atomic number " << AtN(a) ; 
          std::cout << " not currently supported."<< std::endl ;
          break ;

        }
      }

  return basis ;

 } ;

  basis_set build_basis ( common& com, Eigen::Ref<Eigen::VectorXd> AtN, Eigen::Ref<Eigen::MatrixXd> coord) {
  int print = com.prt() ;
  std::string bas_name = com.bnam() ;
  /* Load a basis set */
  const std::string b_cc_pvdz = "cc-pvdz" ;
  const std::string b_6_31g = "6-31g" ;
  const std::string b_sto3g = "sto-3g" ;
  const std::string b_sto1g = "sto-1g" ;
  basis_set basis ;

  if ( bas_name == b_sto3g ){
    basis = load_sto3g( AtN, coord) ;
  } else if ( bas_name == b_sto1g ){
    basis = load_sto1g( AtN, coord) ;
  } else if ( bas_name == b_cc_pvdz ){
    basis = load_cc_pVDZ( AtN, coord) ;
  } else if ( bas_name == b_6_31g ){
    basis = load_6_31g( AtN, coord) ;
    }

  norm_basis ( AtN.size(), basis) ;

  if ( print > 0 ){
    print_basis ( AtN.size(), basis) ;
    }

  return basis ;

    } ;

  void norm_basis ( int natm, basis_set& b) {
  /* Normalize the basis */
    int nsto ;
    double s = 0.0 ;

    for( int i = 0; i < natm; i++ ){
      nsto = b.b[i].nshl ;
      for( int j = 0; j < nsto; j++){
        s = overlap_sto( b.b[i].s[j], b.b[i].c, b.b[i].s[j], b.b[i].c) ;
        b.b[i].s[j].norm = d1/sqrt(s) ;
	}
      } ;

    } ;

  double norm_gprm ( double a, Eigen::Vector3i L) {
  /* Normalize a gaussian primitive */
    double o = d0 ;
    double alpexp = d0 ;
    double s = d1 ;
    double t3_f4 = 3.0e0/4.0e0 ;
 
    for (int i=0; i < 3; i++) {
      s = s*factfact(2*L(i)-1) ;
      alpexp += static_cast<double>(L(i)) ;
      }

    o = d1/sqrt(s) ;
    s = o*pow( 4*a, alpexp/d2) ;
    o = s*pow( (d2*a)/pi, t3_f4) ;

    return o ;

    } ;

  void print_basis( int natm, basis_set& basis){
    /* Dump basis information for debugging */
    int nprm ;
    int nsto ;

    for( int i = 0; i < natm; i++ ){
      std::cout << "Atom number " << i << std::endl ;
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

