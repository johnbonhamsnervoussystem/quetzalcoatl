#include "constants.h"
#include <Eigen/Dense>
#include <iostream>

 template<typename vd, typename vi, class intpre>
  intpre gauprm_ovl( const intpre& za, const intpre& ca, vd&& a, vi&& la, const intpre& zb, const intpre& cb, vd&& b, vi&& lb) {
/*  Two-center overlap integrals evaluated using the Obara-Saika recursion scheme 
    za, zb - Exponents 
    ca, cb - Normalization constants            
    la, lb - Angular momentum information       
    a,  b  - Centers                            
  */
    double thresh = 1.0e-13 ;
    double d1_5 = 3.0e0/2.0e0 ;
    double s, o ;
    intpre zeta, dnrm, xi, rab2 ;
    Eigen::MatrixXd otmp ;
    Eigen::Vector3d p ;
    Eigen::Vector3d t ;
    int na, nb, ia, ib, i ;

    std::cout << " a " << std::endl ;
    zeta = za + zb ;
    xi   = za*zb/zeta ;
    p    = (za*a + zb*b)/zeta ;
    t    = a - b ;
    rab2 = t.dot(t) ;
    std::cout << " rab2 " << std::endl ;
    std::cout << rab2 << std::endl ;
    na   = 1 ;
    nb   = 1 ;
/*    na   = la.maxCoeff() ;
    nb   = lb.maxCoeff() ; */

 /*  Now, we initialize the <s|s> integral */   
      o = pow(pi/zeta, d1_5) * exp(-xi*rab2) ;
      std::cout << " o " << std::endl ;
      std::cout << o << std::endl ;

      std::cout << " b " << std::endl ;
      if(o < thresh) {
        return d0 ;
        }

      std::cout << " d " << std::endl ;
      otmp.resize(na+1,nb+1) ;

/* Increment over each Cartesian */
      std::cout << " e " << std::endl ;
      for (int i = 1; i < 4; i++){
       otmp.setZero() ;
       otmp(1,1) = o ;
       
       std::cout << " g " << std::endl ;
       for (int ib = 0; ib < lb(i-1); ib++){
         otmp(1,ib+2) = (p(i) - b(i))*otmp(1,ib+1) /
          + static_cast<double>(ib)/(d2*zeta)*otmp(1,ib) ;
         }

       std::cout << " h " << std::endl ;
       for (int ib = 0; ib < lb(i-1); ib++){
         for (int ia = 0; ia < la(i-1); ia++){
           otmp(ia+2,ib+1) = (p(i) - a(i))*otmp(ia+1,ib+1) /
            + static_cast<double>(ib)/(d2*zeta)*otmp(ia+1,ib) /
            + static_cast<double>(ia)/(d2*zeta)*otmp(ia,ib+1) ;
           }
         }
             
        std::cout << " f " << std::endl ;
        o = otmp(la(i-1)+1,lb(i-1)+1) ;
        }

      std::cout << " c " << std::endl ;
      otmp.resize( 0, 0) ;

      s = o*ca*cb ;

      return s ;

   } 

int main ( void ){
  /* Testing integral package */
  // Coordinate matrix
  double s ;
  double za=1.0 ;
  double zb=1.0 ;
  double ca=1.0 ;
  double cb=1.0 ;
  
  Eigen::MatrixXd c ;
  Eigen::MatrixXi l ;
  
  l.resize(2,3) ;
  l.setZero() ;

  c.resize(5,3) ;
  c.row(0) << 0.01, 0.00, 0.00 ;
  c.row(1) << 1.01, 0.00, 0.00 ;
  c.row(2) << 2.01, 2.01, 1.01 ;
  c.row(3) << 4.01, 6.01, 4.01 ;
  c.row(4) << -4.01, -6.01, -0.01 ;

  std::cout << " Going into gauprm_ovl " << std::endl ;
  s = gauprm_ovl( za, ca, c.row(0), l.row(0), zb, cb, c.row(1), l.row(1)) ;
  std::cout << s << std::endl ;

  c.resize(0,0) ;
  l.resize(0,0) ;

  return 0 ;

}

