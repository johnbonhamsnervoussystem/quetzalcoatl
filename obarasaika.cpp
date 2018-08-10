#include "basis.h"
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <iostream>
#include <math>
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

  double gauprm_V( double za, double ca, Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3i> la, double zb, double cb, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3i> lb, double q, Eigen::Ref<Eigen::Vector3d> c ) {
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
    Eigen::MatrixXd vmp ;
    Eigen::Tensor< double, 3> vtmp ;
    Eigen::Vector3i lt ;
    Eigen::Vector3d p, abd ;
    int na, nb, ia, ib, i ;
    int ltot ;

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
    o = pow( pi/zeta, d1_5) * exp( -xi*rab2) ;
    v = d0 ;

    vmp.resize( ltot) ;
    vtmp.resize( na+2, nb+2, ltot+1) ;

    abd  = p - c ;
    u = zeta*abd.dot(abd) ;
    for (int m=0; m <= ltot; m++ ) {
      vmp(m) = d2*sqrt( zeta/pi)*o*fboys( m, u) ;
      }

    if ( vmp.cwiseAbs().maxCoeff() < thresh) {
      vmp.resize( 0) ;
      vtmp.resize( 0, 0, 0) ;
      return d0 ;
      }

    for ( int i = 1; i < 4; i++){
      vtmp.setZero() ;
      for( int xx=0; xx =< ltot; xx++){
        vtmp( 1, 1, xx) = vmp( xx) ;
        }

      for ( int ib = 0; ib < lb[i-1]; ib++) {
        for ( m = 0; m < ltot; m++) {
          vtmp( 1, ib+2, m) =  (p(i-1) - b(i-1))*vtmp( 1, ib+1, m) - (p(i-1) - c(i-1))*vtmp( 1, ib+1, m+1)
          + vtmp( 1, ib, m)*static_cast<double>(ib)/(d2*zeta) - vtmp( 1, ib, m+1) ;
          }
        }

      for ( int ib = 0; ib <= lb(i-1); ib++) {
        for ( int ia = 0; ia < la[i-1]; ia++) {
          for ( m = 0; m < ltot; m++) {
            vtmp(ia+2,ib+1,m) = (p(i-1) - a(i-1))*vtmp(ia+1,ib+1,m) - (p(i-1) - c(i-1))*vtmp( ia+1, ib+1, m+1)
            + static_cast<double>(ia)/(d2*zeta)*(vtmp(ia,ib+1,m) - vtmp(ia,ib+1,m+1)) 
            + static_cast<double>(ib)/(d2*zeta)*(vtmp(ia+1,ib,m) - vtmp(ia+1,ib,m+1)) ;
             }
           }
         }

       v = vtmp(la(i),lb(i),:)
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
    ovl.resize( 0, 0) ;

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

    } ;

  double fboys( int m, double t) {
      int k ;
      double gammp, a, f ;
      double r, x, s, fac ;
      const double thresh = 1.0e-4 ;

      a = m + d1/d2 ;
      if ( t < d0 || a < d0) { 
        std::cout << " Negative values in fboys " << std::endl ;
        exit(0) ;
        }

! For very small T, there's a problem with GSer
      if ( t < thresh) {
        r = sqrt(t) ;
        x = d2*t ;
        f = d1/( d2*static_cast<double>(m) + d1) ;
        for ( k = 1,6) {
          f = f + (pow( -t, k)/fact(k))*(d1/( d2*static_cast<double>(m) + d2*static_cast<double>(k) + d1)) ;
          }
      } else if ( t < a + d1) {
        call gser( gammp, a, t) 
        f = gammp/(2*sqrt(t)*t**m)
      } else {
        call gcf( gammp, a, t) ; 
        f = gammp/( d2*sqrt( t)*pow( t, m) ;
      } ;

      return f ;

  } ;

  void gser( double& gammp, double a, double t) {
    const int maxiter = 200 ;
    double const tol = 1.0e-12 ;
    double ap, term, sum ;
    int n ;

    ap   = a ;
    term = d1/ap ;
    sum  = term ;
    for ( n = 1; n <= maxiter; n++){
      ap = ap + d1 ;
      term = term*t/ap ;
      sum = sum + term
      if ( abs(term) < abs(sum)*tol ){
        break ;
        }
      }

    if( abs(term) > abs(sum)*tol ){
      std::cout << " gser needs more terms " << std::endl ;
      call exit(0) ;
      } 

    gammp = sum*exp(-t+a*log(t))

    return ;

  }


  void gcf( double& gammp, double a, double t) {
    const int maxiter = 200 ;
    int m, i ;
    double const tol = 1.0e-12 ;
    double const fpmin = 1.0e-30 ;
    double const d0_5 = d1/d2 ;
    double b, c, d, h, an ;
    double term, gama, dfac ;

    b = t + d1 - a ;
    c = d1/fpmin ;
    d = d1/b ;
    h = d ;

    for ( i = 1; i <= maxiter; i++) {
      an = -static_cast<double>(i)*( static_cast<double>(i) - a) ;
      b  = b + d2 ;
      d  = an*d + b ;
      if ( abs(d) < fpmin) {
        d = fpmin ;
        }
      c  = b + an/c ;
      if ( abs(c) < fpmin) {
        c = fpmin ;
        }
      d = d1/d ;
      term = d*c ;
      h = h*term ;
      if ( abs( term - d1) < tol){
        break ;
        }
    }

    if ( abs( term - d1) > tol) {
      std::cout << " No convergence in gcf " << std::endl ;
      exit(0) ;
     }

    gammp = exp(-t + a*log( t))*h ;

    m = static_cast<int>(round( a - d0_5)) ;
    gama = sqrt(pi)*fact(m)/(pow(2,m)) ;

    gammp = gama - gammp ;

    return ;

  } ;

