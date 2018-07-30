#include "basis.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "obarasaika.h"
#include <vector>

  double gauprm_ovl( const double& za, const double& ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, const double& zb, const double& cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb) {
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

      otmp.resize(na+1,nb+1) ;

/* Increment over each Cartesian */
      for (int i = 1; i < 4; i++){
       otmp.setZero() ;
       otmp(1,1) = o ;
       
       for (int ib = 0; ib < lb[i-1]; ib++){
         otmp(1,ib+2) = (p(i) - b(i))*otmp(1,ib+1) /
          + static_cast<double>(ib)/(d2*zeta)*otmp(1,ib) ;
         }

       for (int ib = 0; ib < lb(i-1); ib++){
         for (int ia = 0; ia < la[i-1]; ia++){
           otmp(ia+2,ib+1) = (p(i) - a(i))*otmp(ia+1,ib+1) /
            + static_cast<double>(ib)/(d2*zeta)*otmp(ia+1,ib) /
            + static_cast<double>(ia)/(d2*zeta)*otmp(ia,ib+1) ;
           }
         }
             
        o = otmp(la[i-1]+1,lb[i-1]+1) ;
        }

      otmp.resize( 0, 0) ;

      s = o*ca*cb ;

      return s ;

   } ;

  double overlap_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb) {
  /* Given two Slater-type Orbitals, return the overlap. */
    int npa = a.nprm ;
    int npb = b.nprm ;
    double s = 0.0 ;

    for ( int i = 0; i < npa; i++) {
      for ( int j = 0; j < npa; j++) {
        s += gauprm_ovl( a.g[i].x, a.g[i].c, ca, a.l, b.g[i].x, b.g[i].c, cb, b.l) ;
        }
      }

    return s ;

    } ;

  void ao_overlap( basis_set& b) {
    /* Given two Slater-type Orbitals, return the overlap. */
    int nbas = b.nbas ;
    int nst1, nst2, jstrt ;
    double s ;
    Eigen::MatrixXd ovl ;
 
    ovl.resize( nbas, nbas) ;
    ovl.setZero() ;
    for ( int i = 0; i < nbas; i++) {
      nst1 = b.b[i].nshl ;
      for ( int j = 0; j < nst1; j++) {
        for ( int k = i; k < nbas; k++) {
          nst2 = b.b[k].nshl ;
          if ( k == i ) {
            jstrt = j ;
          } else {
            jstrt = 0 ;
          }
          for ( int l = jstrt; l < nst2; l++) {
            s = overlap_sto( b.b[i].s[j] , b.b[i].c, b.b[k].s[l] , b.b[k].c) ;
            ovl(k, i) =  s*(b.b[i].s[j].norm)*(b.b[k].s[l].norm) ;
            ovl(i, k) = ovl(k, i) ;
            }
          }
        }
      }

    std::cout << ovl << std::endl ;

    } ;

