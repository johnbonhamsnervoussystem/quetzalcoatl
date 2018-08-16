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
      std::cout << " m " << m << std::endl ;
      std::cout << " u " << u << std::endl ;
      std::cout << "fboys( m, u) = " << fboys( m, u) << std::endl ;
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
      Real (Kind=pr)              :: ERInt
! Local variables
      Real (Kind=pr), Allocatable :: ETmp(:,:,:,:,:), E(:)
      Real (Kind=pr) :: Zeta, Eta, Rho, DNrm, Oab, Ocd
      Eigen::Vector3d p, q, w ;

/*      Real (Kind=pr) :: P(3), Q(3), W(3)
      Real (Kind=pr) :: Rab2, Rcd2, Rpq2, T
      Integer :: NA, NB, NC, ND, LTot
      Integer :: IA, IB, IC, ID, M, I */
   
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

      zeta = za + zb 
      eta  = zc + zd
      rho  = eta*zeta/(eta+zeta)
      P    = (ZA*A + ZB*B)/Zeta
      Q    = (ZC*C + ZD*D)/Eta
      W    = (Zeta*P + Eta*Q)/(Zeta+Eta)
      Rab2 = Dot_Product(A-B,A-B)
      Rcd2 = Dot_Product(C-D,C-D)
      Rpq2 = Dot_Product(P-Q,P-Q)
      T    = Rho*Rpq2
      LTot = Sum(LA+LB+LC+LD)
      NA   = MaxVal(LA)
      NB   = MaxVal(LB)
      NC   = MaxVal(LC)
      ND   = MaxVal(LD)


!=====================================================!
!  After initializing the <s|s> integrals, we start   !
!  incrementing the Cartesian components of each      !
!  function in turn.  Since increasing the angular    !
!  momentum requires auxiliary integrals indexed by   !
!  m, we have to loop over m as well.                 !
!=====================================================!

      Oab = (Pi/Zeta)**F32 * Exp(-ZA*ZB*Rab2/Zeta)
      Ocd = (Pi/ Eta)**F32 * Exp(-ZC*ZD*Rcd2/Eta)
      Allocate(E(0:LTot), ETmp(-1:NA,-1:NB,-1:NC,-1:ND,0:LTot))

      Do M = 0,LTot
       E(M) = Two*Sqrt(Rho/Pi)*Oab*Ocd*FBoys(M,T)
      End Do
! TMH: New!
      If(MaxVal(Abs(E(0:LTot))) < Thresh) Then
        ERInt = Zero
        Deallocate(E, ETmp)
        Return
      End If
! TMH: Done!

      Do I = 1,3
       ETmp = Zero
       ETmp(0,0,0,0,:) = E

       Do ID = 0,LD(I)-1
        Do M = 0,LTot-1
         ETmp(0,0,0,ID+1,M) = (Q(I) - D(I))*ETmp(0,0,0,ID,M)             &
                            + (W(I) - Q(I))*ETmp(0,0,0,ID,M+1)           &
               + Real(ID)/(Two*Eta)*(ETmp(0,0,0,ID-1,M)                  &
                                   - ETmp(0,0,0,ID-1,M+1)*Rho/Eta)
        End Do
       End Do

       Do ID = 0,LD(I)
        Do IC = 0,LC(I)-1
         Do M = 0,LTot-1
          ETmp(0,0,IC+1,ID,M) = (Q(I) - C(I))*ETmp(0,0,IC,ID,M)          &
                              + (W(I) - Q(I))*ETmp(0,0,IC,ID,M+1)        &
               + Real(ID)/(Two*Eta)*(ETmp(0,0,IC,ID-1,M)                 &
                                   - ETmp(0,0,IC,ID-1,M+1)*Rho/Eta)      &
               + Real(IC)/(Two*Eta)*(ETmp(0,0,IC-1,ID,M)                 &
                                   - ETmp(0,0,IC-1,ID,M+1)*Rho/Eta)
         End Do
        End Do
       End Do

       Do ID = 0,LD(I)
        Do IC = 0,LC(I)
         Do IB = 0,LB(I)-1
          Do M = 0,LTot-1
           ETmp(0,IB+1,IC,ID,M) = (P(I) - B(I))*ETmp(0,IB,IC,ID,M)       &
                                + (W(I) - P(I))*ETmp(0,IB,IC,ID,M+1)     &
               + Real(ID)/(Two*Eta+Two*Zeta)*ETmp(0,IB,IC,ID-1,M+1)      &
               + Real(IC)/(Two*Eta+Two*Zeta)*ETmp(0,IB,IC-1,ID,M+1)      &
               + Real(IB)/(Two*Zeta)*(ETmp(0,IB-1,IC,ID,M)               &
                                    - ETmp(0,IB-1,IC,ID,M+1)*Rho/Zeta)
          End Do
         End Do
        End Do
       End Do

       Do ID = 0,LD(I)
        Do IC = 0,LC(I)
         Do IB = 0,LB(I)
          Do IA = 0,LA(I)-1
           Do M = 0,LTot-1
            ETmp(IA+1,IB,IC,ID,M) = (P(I) - A(I))*ETmp(IA,IB,IC,ID,M)    &
                                  + (W(I) - P(I))*ETmp(IA,IB,IC,ID,M+1)  &
               + Real(ID)/(Two*Eta+Two*Zeta)*ETmp(IA,IB,IC,ID-1,M+1)     &
               + Real(IC)/(Two*Eta+Two*Zeta)*ETmp(IA,IB,IC-1,ID,M+1)     &
               + Real(IB)/(Two*Zeta)*(ETmp(IA,IB-1,IC,ID,M)              &
                                    - ETmp(IA,IB-1,IC,ID,M+1)*Rho/Zeta)  &
               + Real(IA)/(Two*Zeta)*(ETmp(IA-1,IB,IC,ID,M)              &
                                    - ETmp(IA-1,IB,IC,ID,M+1)*Rho/Zeta)
           End Do
          End Do
         End Do
        End Do
       End Do
       E = ETmp(LA(I),LB(I),LC(I),LD(I),:)
      End Do

      ERInt = E(0)*DNrm
      DeAllocate(E,ETmp)

      Return
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

    } ;

  double fboys( int k, double t) {
      double f, m, n ;
      /*  Let's not mess around.  Just compute this*/

      m = static_cast<double>(k) + d1/d2 ;
      n = static_cast<double>(k) + 3.0e0/d2 ;
      f = gsl_sf_hyperg_1F1( m, n, -t)/( d2*static_cast<double>(k) + d1) ;
      
      return f ;

    } ;

