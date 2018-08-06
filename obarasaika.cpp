#include "basis.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "obarasaika.h"
#include <vector>

  double gauprm_ovl( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) {
/*  Two-center overlap integrals evaluated using the Obara-Saika recursion scheme 
    za, zb - Exponents 
    ca, cb - Normalization constants            
    la, lb - Angular momentum information       
    a,  b  - Centers                            
  */

    double thresh = 1.0e-13 ;
    double d1_5 = 3.0e0/2.0e0 ;
    double s, o ;
    double zeta, dnrm, xi, rab2 ;
    Eigen::MatrixXd otmp ;
    Eigen::Vector3d p ;
    Eigen::Vector3d t ;
    int na, nb, ia, ib, i ;

    zeta = za + zb ;
    xi   = za*zb/zeta ;
    p    = (za*a + zb*b)/zeta ;
    t    = a - b ;
    rab2 = t.dot(t) ;
    na   = la.maxCoeff() ;
    nb   = lb.maxCoeff() ; 

 /*  Now, we initialize the <s|s> integral */   
    o = pow(pi/zeta, d1_5) * exp(-xi*rab2) ;

    if(o < thresh) {
      return d0 ;
    } else if ( na == 0 && nb == 0 ) {
      s = o*ca*cb ;
      return s ;
      }

    otmp.resize(na+2,nb+2) ;

/* Increment over each Cartesian */
    for (int i = 1; i < 4; i++){
      otmp.setZero() ;
      otmp(1,1) = o ;
      for ( int ib = 0; ib < lb[i-1]; ib++){
        otmp(1,ib+2) = (p(i-1) - b(i-1))*otmp(1,ib+1) 
         + otmp(1,ib)*static_cast<double>(ib)/(d2*zeta);
        }

      for ( int ib = 0; ib < lb(i-1)+1; ib++){
        for ( int ia = 0; ia < la[i-1]; ia++){
          otmp(ia+2,ib+1) = (p(i-1) - a(i-1))*otmp(ia+1,ib+1) 
           + static_cast<double>(ib)/(d2*zeta)*otmp(ia+1,ib) 
           + static_cast<double>(ia)/(d2*zeta)*otmp(ia,ib+1) ;
          }
        }
 
       o = otmp(la[i-1]+1,lb[i-1]+1) ;
       }

      otmp.resize( 0, 0) ;

      s = o*ca*cb ;

      return s ;

    } ;


  double gauprm_T( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) {
/*  Two-center kinetic energy integrals evaluated using the Obara-Saika recursion scheme 
    za, zb - Exponents 
    ca, cb - Normalization constants            
    la, lb - Angular momentum information       
    a,  b  - Centers                            
  */

    double thresh = 1.0e-13 ;
    double d1_5 = 3.0e0/2.0e0 ;
    double k, o, t ;
    double zeta, dnrm, xi, rab2 ;
    Eigen::MatrixXd otmp, ttmp ;
    Eigen::Vector3d p, abd ;
    int na, nb, ia, ib, i ;

    zeta = za + zb ;
    xi   = za*zb/zeta ;
    p    = (za*a + zb*b)/zeta ;
    abd    = a - b ;
    rab2 = abd.dot(abd) ;
    na   = la.maxCoeff() ;
    nb   = lb.maxCoeff() ; 

 /*  Now, we initialize the <s|s> integral */   
    o = pow(pi/zeta, d1_5) * exp(-xi*rab2) ;
    t = xi*(3.0e0 - 2.0e0*xi*rab2)*o ;

    if ( o < thresh) {
      return d0 ;
      }

    otmp.resize(na+2,nb+2) ;
    ttmp.resize(na+2,nb+2) ;

    for (int i = 1; i < 4; i++){
      o = d0 ;
      t = d0 ;
      otmp(1,1) = o ;
      ttmp(1,1) = t ;

      for ( int ib = 0; ib < lb[i-1]; ib++){
        otmp(1,ib+2) = (p(i-1) - b(i-1))*otmp(1,ib+1) 
         + otmp(1,ib)*static_cast<double>(ib)/(d2*zeta);
        ttmp(1,ib+2) = (p(i-1) - b(i-1))*ttmp(1,ib+1) + d2*xi*otmp(1,ib+2)
          + static_cast<double>(ib)*(ttmp(1,ib)/(d2*zeta) - otmp(1,ib)*xi/zb) ;
         }

      for ( int ib = 0; ib < lb(i-1)+1; ib++){
        for ( int ia = 0; ia < la[i-1]; ia++){
         otmp(ia+2,ib+1) = (p(i-1) - a(i-1))*otmp(ia+1,ib+1)
           + otmp(ia+1,ib)*static_cast<double>(ib)/(d2*zeta) 
           + otmp(ia,ib+1)*static_cast<double>(ia)/(d2*zeta) ;
         ttmp(ia+2,ib+1) = (p(i-1) - a(i-1))*ttmp(ia+1,ib+1)   
                       + d2*xi*otmp(ia+2,ib+1)        
           + static_cast<double>(ia)*(ttmp(ia,ib+1)/(d2*zeta) 
                                 - otmp(ia,ib+1)*xi/za) 
           + static_cast<double>(ib)*ttmp(ia+1,ib)/(d2*zeta) ;
          }
        }

       o = otmp(la[i-1]+1,lb[i-1]+1) ;
       t = ttmp(la[i-1]+1,lb[i-1]+1) ;
       }

      otmp.resize( 0, 0) ;
      ttmp.resize( 0, 0) ;

      k = t*ca*cb ;
 
      return k ;

  } ;


  double overlap_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) {
  /* Given two Slater-type Orbitals, return the overlap. */
    int npa = a.nprm ;
    int npb = b.nprm ;
    double s = 0.0 ;

    for ( int i = 0; i < npa; i++) {
      for ( int j = 0; j < npb; j++) {
        s += gauprm_ovl( a.g[i].x, a.g[i].c, ca, a.l, b.g[j].x, b.g[j].c, cb, b.l) ;
        }
      }

    return s ;

    } ;

  double kinetic_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) {
  /* Given two Slater-type Orbitals, return the kinetic energy. */
    int npa = a.nprm ;
    int npb = b.nprm ;
    double t = 0.0 ;

    for ( int i = 0; i < npa; i++) {
      for ( int j = 0; j < npb; j++) {
        t += gauprm_T( a.g[i].x, a.g[i].c, ca, a.l, b.g[j].x, b.g[j].c, cb, b.l) ;
        }
      }

    return t ;

    } ;

  void ao_overlap( int natm, basis_set& b) {
    /* Given a Slater-type Orbital basis, return the overlap. */
    int nbas = b.nbas ;
    int ind = -1 ;
    int jnd ;
    int nst1, nst2, jstrt ;
    double s ;
    Eigen::MatrixXd ovl ;
 
    ovl.resize( nbas, nbas) ;
    ovl.setZero() ;
    for ( int i = 0; i < natm; i++) {
      nst1 = b.b[i].nshl ;
      for ( int j = 0; j < nst1; j++) {
        ind ++ ;
        jnd = ind - 1 ;
        for ( int k = i; k < natm; k++) {
          nst2 = b.b[k].nshl ;
          if ( k == i ) {
            jstrt = j ;
          } else {
            jstrt = 0 ;
          }
          for ( int l = jstrt; l < nst2; l++) {
            jnd ++ ;
            s = overlap_sto( b.b[i].s[j] , b.b[i].c, b.b[k].s[l] , b.b[k].c) ;
            ovl( ind, jnd) =  s*(b.b[i].s[j].norm)*(b.b[k].s[l].norm) ;
            ovl( jnd, ind) = ovl( ind, jnd) ;
            }
          }
        }
      }

    std::cout << ovl << std::endl ;

    } ;

  void ao_kinetic( int natm, basis_set& b) {
    /* Given a Slater-type Orbital basis, return the kinetic
       energy. */
    int nbas = b.nbas ;
    int ind = -1 ;
    int jnd ;
    int nst1, nst2, jstrt ;
    double k ;
    Eigen::MatrixXd T ;
 
    T.resize( nbas, nbas) ;
    T.setZero() ;
    for ( int i = 0; i < natm; i++) {
      nst1 = b.b[i].nshl ;
      for ( int j = 0; j < nst1; j++) {
        ind ++ ;
        jnd = ind - 1 ;
        for ( int k = i; k < natm; k++) {
          nst2 = b.b[k].nshl ;
          if ( k == i ) {
            jstrt = j ;
          } else {
            jstrt = 0 ;
          }
          for ( int l = jstrt; l < nst2; l++) {
            jnd ++ ;
            k = kinetic_sto( b.b[i].s[j] , b.b[i].c, b.b[k].s[l] , b.b[k].c) ;
            T( ind, jnd) =  k ;
            T( jnd, ind) = T( ind, jnd) ;
            }
          }
        }
      }

    std::cout << T << std::endl ;

    } ;

