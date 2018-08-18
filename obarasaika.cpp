#include "basis.h"
#include <cmath>
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <gsl/gsl_sf_hyperg.h>
#include <iostream>
#include "obarasaika.h"
#include "util.h"
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
    double zeta, xi, rab2 ;
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

 /* <s|s> */   
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
    double zeta, xi, rab2 ;
    Eigen::MatrixXd otmp, ttmp ;
    Eigen::Vector3d p, abd ;
    int na, nb, ia, ib, i ;

    zeta = za + zb ;
    xi   = za*zb/zeta ;
    p    = (za*a + zb*b)/zeta ;
    abd  = a - b ;
    rab2 = abd.dot(abd) ;
    na   = la.maxCoeff() ;
    nb   = lb.maxCoeff() ; 

/* <s|s> */
    o = pow(pi/zeta, d1_5) * exp(-xi*rab2) ;
    t = xi*(3.0e0 - 2.0e0*xi*rab2)*o ;

    if ( o < thresh) {
      return d0 ;
      }

    otmp.resize(na+2,nb+2) ;
    ttmp.resize(na+2,nb+2) ;

    for (int i = 1; i < 4; i++){
      otmp.setZero() ;
      ttmp.setZero() ;
      otmp(1,1) = o ;
      ttmp(1,1) = t ;

      for ( int ib = 0; ib < lb[i-1]; ib++) {
        otmp(1,ib+2) = (p(i-1) - b(i-1))*otmp(1,ib+1) 
         + otmp(1,ib)*static_cast<double>(ib)/(d2*zeta);
        ttmp(1,ib+2) = (p(i-1) - b(i-1))*ttmp(1,ib+1) + d2*xi*otmp(1,ib+2)
          + static_cast<double>(ib)*(ttmp(1,ib)/(d2*zeta) - otmp(1,ib)*xi/zb) ;
         }

      for ( int ib = 0; ib < lb(i-1)+1; ib++) {
        for ( int ia = 0; ia < la[i-1]; ia++) {
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

  double gauprm_V( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb, Eigen::Ref<Eigen::Vector3d> c, double q ) { 
/*  Two-center potential energy integrals evaluated using the Obara-Saika recursion scheme 
    za, zb - Exponents 
    ca, cb - Normalization constants            
    la, lb - Angular momentum information       
    a,  b  - Centers                            
    q      - charge
    c      - Nuclear Coordinates
  */

    double thresh = 1.0e-13 ;
    double d1_5 = 3.0e0/2.0e0 ;
    double v, o, u ;
    double zeta, xi, rab2 ;
    Eigen::VectorXd vmp ;
    Eigen::Tensor< double, 3> vtmp( 1, 1, 1) ;
    Eigen::Vector3i lt ;
    Eigen::Vector3d p, abd ;
    int na, nb, ia, ib, i ;
    int ltot, m ;

    zeta = za + zb ;
    xi   = za*zb/zeta ;
    p    = (za*a + zb*b)/zeta ;
    abd  = a - b ;
    rab2 = abd.dot(abd) ;
    lt   = la + lb ;
    ltot = lt.sum() ;
    na   = la.maxCoeff() ;
    nb   = lb.maxCoeff() ; 

/* <s|s> */
    o = pow( pi/zeta, d1_5)*exp( -xi*rab2) ;
    v = d0 ;

    vmp.resize( ltot+1) ;
    vtmp = Eigen::Tensor <double, 3>( na+2, nb+2, ltot+1) ;

    abd = p - c ;
    u = zeta*( abd.dot( abd)) ;
    for ( int m = 0; m <= ltot; m++ ) {
      vmp(m) = d2*sqrt( zeta/pi)*o*fboys( m, u) ;
      }

    if ( vmp.cwiseAbs().maxCoeff() < thresh) {
      vmp.resize( 0) ;
      vtmp = Eigen::Tensor <double, 3>( 0, 0, 0) ;
      return d0 ;
      }

    for ( int i = 1; i < 4; i++){
      vtmp.setZero() ;
      for( int xx=0; xx <= ltot; xx++){
        vtmp( 1, 1, xx) = vmp( xx) ;
        }

      for ( int ib = 0; ib < lb[i-1]; ib++) {
        for ( m = 0; m < ltot; m++) {
          vtmp( 1, ib+2, m) = (p(i-1) - b(i-1))*vtmp( 1, ib+1, m) - (p(i-1) - c(i-1))*vtmp( 1, ib+1, m+1)
          + (vtmp( 1, ib, m) - vtmp( 1, ib, m+1))*static_cast<double>(ib)/(d2*zeta) ;
          }
        }

      for ( int ib = 0; ib <= lb(i-1); ib++) {
        for ( int ia = 0; ia < la[i-1]; ia++) {
          for ( m = 0; m < ltot; m++) {
            vtmp( ia+2, ib+1, m) = (p(i-1) - a(i-1))*vtmp( ia+1, ib+1, m) - (p(i-1) - c(i-1))*vtmp( ia+1, ib+1, m+1)
            + (vtmp(ia,ib+1,m) - vtmp(ia,ib+1,m+1))*static_cast<double>(ia)/(d2*zeta) 
            + (vtmp(ia+1,ib,m) - vtmp(ia+1,ib,m+1))*static_cast<double>(ib)/(d2*zeta) ;
             }
           }
         }

      for( int xx=0; xx <= ltot; xx++){
        vmp( xx) = vtmp( la[i-1]+1, lb[i-1]+1, xx) ;
        }
       }

      v = vmp(0)*ca*cb*q ;

      vmp.resize( 0) ;
      vtmp = Eigen::Tensor <double, 3>( 0, 0, 0) ;
 
      return v ;

  } ;

  double gauprm_r12( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb, double zc, double cc, Eigen::Ref<Eigen::Vector3d> c, Eigen::Ref<Eigen::Vector3i> lc, double zd, double cd, Eigen::Ref<Eigen::Vector3d> d, Eigen::Ref<Eigen::Vector3i> ld )  {
/*  Two electron integrals evaluated using the Obara-Saika recursion scheme 
    za, zb, zc, zd - Exponents 
    ca, cb, cc, cd - Normalization constants            
    la, lb, lc, ld - Angular momentum information       
    a, b, c, d - Centers                            
  */
    int na, nb, nc, nd, ltot ;
    int ia, ib, ic, id, m, i, xx ;
    double thresh = 1.0e-13 ;
    double d1_5 = 3.0e0/2.0e0 ;
    double rab2, rcd2, rpq2, t, er12 ;
    double zeta, eta, rho, oab, ocd ;
    Eigen::Vector3i lt ;
    Eigen::Vector3d p, q, w, abd ;
    Eigen::VectorXd e ;
    Eigen::Tensor< double, 5> etmp( 1, 1, 1, 1, 1) ;

/*   
!====================================================!
!  Electron repulsion integrals evaluated using the  !
!  Obara-Saika recursion scheme.                     !    
!                                                    !
!  The functions we're integrating have:             !
!    ZK = Exponent of function K                     !
!    CK = Normalization constant of function K       !
!    LK = Angular momentum for function K            !
!    K  = Vector pointing to center K                !
!  Note that this integral is in Mulliken notation,  !
!  that is, with functions A, B, C, and D we have    !
!    ERI(A,B,C,D) = (A B | 1/r12 | C D)              !
!                 = <A C | 1/r12 | B D>.             !
!====================================================!
*/
    zeta = za + zb ;
    eta  = zc + zd ;
    rho  = eta*zeta/(eta+zeta) ;
    p    = (za*a + zb*b)/zeta ;
    q    = (zc*c + zd*d)/eta ;
    w    = (zeta*p + eta*q)/(zeta+eta) ;
    abd = a - b ;
    rab2 = abd.dot(abd) ;
    abd = c - d ;
    rcd2 = abd.dot(abd) ;
    abd = p - q ;
    rpq2 = abd.dot(abd) ;
    t = rho*rpq2 ;
    lt = la + lb + lc + ld ;
    ltot = lt.sum() ;
    na = la.maxCoeff() ;
    nb = lb.maxCoeff() ; 
    nc = lc.maxCoeff() ;
    nd = ld.maxCoeff() ; 
/*
!=====================================================!
!  After initializing the <s|s> integrals, we start   !
!  incrementing the Cartesian components of each      !
!  function in turn.  Since increasing the angular    !
!  momentum requires auxiliary integrals indexed by   !
!  m, we have to loop over m as well.                 !
!=====================================================!
*/
    oab = pow( pi/zeta, d1_5)*exp( -za*zb*rab2/zeta) ;
    ocd = pow( pi/eta, d1_5)*exp( -zc*zd*rcd2/eta) ;

    e.resize( ltot+1) ;
    etmp = Eigen::Tensor <double, 5>( na+2, na+2, na+2, na+2, ltot+1) ;

    for( m = 0; m <= ltot; m++) {
      e(m) = d2*sqrt( rho/pi)*oab*ocd*fboys( m, t) ;
      }

    if ( e.cwiseAbs().maxCoeff() < thresh) {
      e.resize( 0) ;
      etmp = Eigen::Tensor <double, 5>( 0, 0, 0, 0, 0) ;
      return d0 ;
      }

    for( i = 1; i < 4; i++) {
      etmp.setZero() ;

      for( xx=0; xx <= ltot; xx++) {
        etmp( 1, 1, 1, 1, xx) = e( xx) ;
        }

       for( id = 0; id < ld( i-1); id++) {
         for( m = 0; m < ltot; m++) {
           etmp( 1, 1, 1, id+2, m) = (q(i-1) - d(i-1))*etmp( 1, 1, 1, id+1, m) 
           + (w(i-1) - q(i-1))*etmp( 1, 1, 1, id+1, m+1)
           + static_cast<double>(id)/(d2*eta)*(etmp( 1, 1, 1, id, m) 
           - etmp( 1, 1, 1, id, m+1)*rho/eta) ;
           }
         }

       for( id = 0; id <= ld( i-1); id++) {
         for( ic = 0; ic < lc( i-1); ic++) {
           for( m = 0; m < ltot; m++) {
             etmp( 1, 1, ic+2, id+1, m) = (q(i-1) - c(i-1))*etmp( 1, 1, ic+1, id+1, m)
               + (w(i-1) - q(i-1))*etmp( 1, 1, ic+1, id+1, m+1)  
               + static_cast<double>(id)/(d2*eta)*(etmp( 1, 1, ic+1, id, m) 
               - etmp( 1 , 1, ic+1, id, m+1)*rho/eta) + static_cast<double>(ic)/(d2*eta)*(etmp( 1, 1, ic, id+1, m)
               - etmp( 1, 1, ic, id+1, m+1)*rho/eta) ;
             }
           }
         }

       for( id = 0; id <= ld( i-1); id++) {
         for( ic = 0; ic <= lc( i-1); ic++) {
           for( ib = 0; ib < lb( i-1); ib++) {
             for( m = 0; m < ltot; m++) {
               etmp( 1, ib+2, ic+1, id+1, m) = (p(i-1) - b(i-1))*etmp( 1, ib+1, ic+1, id+1, m)
               + (w(i-1) - p(i-1))*etmp( 1,ib+1, ic+1, id+1, m+1) 
               + etmp( 1, ib+1, ic+1, id, m+1)*static_cast<double>(id)/( d2*eta + d2*zeta)
               + etmp( 1, ib+1, ic, id+1, m+1)*static_cast<double>(ic)/( d2*eta + d2*zeta)
               + static_cast<double>(ib)/(d2*zeta)*(etmp( 1, ib, ic+1, id+1, m)
               - etmp( 1, ib, ic+1, id+1, m+1)*rho/zeta) ;
               }
             }
           }
         }

       for( id = 0; id <= ld( i-1); id++) {
         for( ic = 0; ic <= lc( i-1); ic++) {
           for( ib = 0; ib <= lb( i-1); ib++) {
             for( ia = 0; ia < la( i-1); ia++) {
               for( m = 0; m < ltot; m++) {
                 etmp( ia+2, ib+1, ic+1, id+1, m) = (p(i-1) - a(i-1))*etmp( ia+1, ib+1, ic+1, id+1, m)
                 + (w(i-1) - p(i-1))*etmp( ia+1, ib+1, ic+1, id+1, m+1)
               + etmp( ia+1, ib+1, ic+1, id, m+1)*static_cast<double>(id)/( d2*eta + d2*zeta)
               + etmp( ia+1, ib+1, ic, id+1, m+1)*static_cast<double>(ic)/( d2*eta + d2*zeta)
               + (etmp( ia+1, ib, ic+1, id+1, m) - etmp( ia+1, ib, ic+1, id+1, m+1)*rho/zeta)*static_cast<double>(ib)/(d2*zeta) 
               + ( etmp( ia, ib+1, ic+1, id+1, m) - etmp( ia, ib+1, ic+1, id+1, m+1)*rho/zeta) *static_cast<double>(ia)/(d2*zeta) ;
                 }
               }
             }
           }
         }

      for( int xx=0; xx <= ltot; xx++) {
        e( xx) = etmp( la(i-1)+1, lb(i-1)+1, lc(i-1)+1, ld(i-1)+1, xx) ;
        }

      }

      er12 = e(0)*ca*cb*cc*cd ;
      e.resize( 0) ;
      etmp = Eigen::Tensor <double, 5>( 0, 0, 0, 0, 0) ;

      return er12 ;

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

  double nucelecV_sto( sto& a, Eigen::Ref<Eigen::Vector3d> ca, sto& b, Eigen::Ref<Eigen::Vector3d> cb, Eigen::Ref<Eigen::MatrixXd> n_c, Eigen::Ref<Eigen::VectorXd> q) {
  /* Given two Slater-type Orbitals, return the kinetic energy. */
    int npa = a.nprm ;
    int npb = b.nprm ;
    double v = 0.0 ;
    Eigen::Vector3d cn ;

    for( int k = 0; k < q.size(); k++) {
      cn = n_c.row(k) ;
      for ( int i = 0; i < npa; i++) {
        for ( int j = 0; j < npb; j++) {
          v += gauprm_V( a.g[i].x, a.g[i].c, ca, a.l, b.g[j].x, b.g[j].c, cb, b.l, cn, q(k)) ;
          }
        }
      }

    return v ;

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
    ovl.resize( 0, 0) ;

    return ;

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
            T( ind, jnd) = kinetic_sto( b.b[i].s[j] , b.b[i].c, b.b[k].s[l] , b.b[k].c) ;
            T( jnd, ind) = T( ind, jnd) ;
            }
          }
        }
      }

    std::cout << T << std::endl ;
    T.resize( 0, 0) ;

    return ;

    } ;

  void ao_eN_V( int natm, basis_set& b, Eigen::Ref<Eigen::MatrixXd> n_c, Eigen::Ref<Eigen::VectorXd> q ) {
    /* Given a Slater-type Orbital basis and the nuclear information, return
    the electron nuclear potential integrals */
    int nbas = b.nbas ;
    int ind = -1 ;
    int jnd ;
    int nst1, nst2, jstrt ;
    double k ;
    Eigen::MatrixXd V ;
 
    V.resize( nbas, nbas) ;
    V.setZero() ;

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
            V( ind, jnd) = nucelecV_sto( b.b[i].s[j] , b.b[i].c, b.b[k].s[l] , b.b[k].c, n_c, q) ;
            V( jnd, ind) = V( ind, jnd) ;
            }
          }
        }
      }

    std::cout << V << std::endl ;
    V.resize( 0, 0) ;

    return ;

    } ;

  void ao_tei( int natm, basis_set& b ) {
    /* Given a Slater-type Orbital basis and the nuclear information, return
    the two electron repulsion integrals over Muliken ao indexes */
    int nbas = b.nbas ;
    int ind = -1 ;
    int jnd ;
    int nst1, nst2, jstrt ;
    double k ;
    Eigen::MatrixXd V ;
 
    V.resize( nbas, nbas) ;
    V.setZero() ;

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
            V( ind, jnd) = nucelecV_sto( b.b[i].s[j] , b.b[i].c, b.b[k].s[l] , b.b[k].c, n_c, q) ;
            V( jnd, ind) = V( ind, jnd) ;
            }
          }
        }
      }

    std::cout << V << std::endl ;
    V.resize( 0, 0) ;

    return ;

    } ;


