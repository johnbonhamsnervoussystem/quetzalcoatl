#include "constants.h"
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "binio.h"
#include "evalm.h"
#include "qtzio.h"
#include "tei.h"

/* 
   evalm is a collection of routines that evaluates matrix elements.
   I think these routines are not properly accounting for the complex
   conjugation of densities.  I will have to look at this at some point.
*/

void ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> p, 
  Eigen::Ref<Eigen::MatrixXd> g, const int nb) {

/* 
  G is split into four quadrants.  

  | alpha alpha | alpha beta |
  ---------------------------- 
  | beta alpha  | beta beta  | 
  
  For restricted Hartree-Fock calculations, the orbitals are either
  alpha or beta, not a combination of both so the mixed blocks are
  zero.  Further, the alpha and beta blocks are identical so we only
  need to calculate the coulomb and exchange for a single block.
  
  G_{mu nu} = sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha + 
  p_{sigma lambda}^{beta beta })( mu nu| lambda sigma) - 
  sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha)( mu nu| sigma lambda ) 
  
*/

    Eigen::MatrixXd t_p ;

    g.setZero() ;
    t_p.resize( nb, nb) ;

    t_p = d2*p ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, t_p, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    t_p.resize( 0, 0) ;

    return ;

} ;

void ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> g, const int nb) {

 /*  
  * This is an overloaded function with complex matrices for complex
  * restricted Hartree-Fock.
  */

    Eigen::MatrixXcd t_p ;

    g.setZero() ;
    t_p.resize( nb, nb) ;

    t_p = d2*p ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, t_p, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    t_p.resize( 0, 0) ;

    return ;
} ;

void ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> pt, 
     Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> g, const int nb) {

 /* G is split into four quadrants.
  *
  *    | alpha alpha | alpha beta |
  *    ----------------------------
  *    | beta alpha  | beta beta  |
  *
  *    For urestricted Hartree-Fock calculations, the orbitals are either
  *    alpha or beta, not a combination of both so the mixed blocks are
  *    zero.  However the density of both sets of MOs are needed to 
  *    calculate the coulomb repulsion while only the same spin is needed
  *    for exchange.
  *
  *    G_{mu nu} = sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha + 
  *    p_{sigma lambda}^{beta beta })( mu nu| lambda sigma) - 
  *    sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha)( mu nu| sigma lambda ) 
  *
  */

    g.setZero() ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, pt, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    return ;
} ;

void ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> pt, 
     Eigen::Ref<Eigen::MatrixXcd> p, Eigen::Ref<Eigen::MatrixXcd> g, const int nb) {

 /*  
  * This is an overloaded function with complex matrices for complex
  * unrestricted Hartree-Fock.
  */
    g.setZero() ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, pt, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    return ;
} ;

void ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> g, const int nb) {

 /* G is split into four quadrants.  
  *
  *    | alpha alpha | alpha beta |
  *    ---------------------------- 
  *    | beta alpha  | beta beta  | 
  *
  *    The ghf orbitals are a combination of alpha and beta spins. In this case
  *    the mixed spin blocks have a contribution too. The same spin blocks have 
  *    the form
  *
  *    G_{mu nu} = sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha + 
  *    p_{sigma lambda}^{beta beta })( mu nu| lambda sigma) - 
  *    sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha)( mu nu| sigma lambda ) 
  *
  *    The mixed spin block doesn't have a coulomb term but it has an exchange 
  *    term.  
  *
  *    G_{mu nu}^{alpha beta) = -sum_{lambda sigma}(p_{sigma lambda}^{alpha beta)*
  *    ( mu nu| sigma lambda) 
  *
  *    Pass each quadrant to be contracted with the two electron integrals.
  *
  */

    Eigen::MatrixXd t_p ;

    g.setZero() ;
    t_p.resize( nb, nb) ;

    t_p.block(0,0,nb,nb) = p.block(0,0,nb,nb) + p.block(nb,nb,nb,nb) ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, t_p, g.block( 0, 0, nb, nb), nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p.block(0,0,nb,nb), g.block( 0, 0, nb, nb), nb) ;

    /* Do the coulomb terms for the beta beta block  */
    coulblk( intarr, t_p, g.block( nb, nb, nb, nb), nb) ;

    /* Do the exchange terms for the beta beta block  */
    exchblk( intarr, p.block(nb,nb,nb,nb), g.block( nb, nb, nb, nb), nb) ;

    /* Do the exchange terms for the alpha beta block  */
    exchblk( intarr, p.block(0,nb,nb,nb), g.block( 0, nb, nb, nb), nb) ;

    /* Do the exchange terms for the beta alpha block  */
    exchblk( intarr, p.block(nb,0,nb,nb), g.block( nb, 0, nb, nb), nb) ;

    t_p.resize( 0, 0) ;

    return ;
} ;

void ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> g, const int nb) {

 /*   
  * This is the complex overloaded version of ctr2eg
  */
    Eigen::MatrixXcd t_p ;

    g.setZero() ;
    t_p.resize( nb, nb) ;

    t_p = p.block(0,0,nb,nb) + p.block(nb,nb,nb,nb) ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, t_p, g.block( 0, 0, nb, nb), nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p.block(0,0,nb,nb), g.block( 0, 0, nb, nb), nb) ;

    /* Do the coulomb terms for the beta beta block  */
    coulblk( intarr, t_p, g.block( nb, nb, nb, nb), nb) ;

    /* Do the exchange terms for the beta beta block  */
    exchblk( intarr, p.block(nb,nb,nb,nb), g.block( nb, nb, nb, nb), nb) ;

    /* Do the exchange terms for the alpha beta block  */
    exchblk( intarr, p.block(0,nb,nb,nb), g.block( 0, nb, nb, nb), nb) ;

    /* Do the exchange terms for the beta alpha block  */
    exchblk( intarr, p.block(nb,0,nb,nb), g.block( nb, 0, nb, nb), nb) ;

    t_p.resize( 0, 0) ;

    return ;
} ;

  void ctrPairg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> k, 
    Eigen::Ref<Eigen::MatrixXd> D, const int nb) {

 /* 
  *    This routine contracts the two electron integrals to build the pairing
  *    field for Hartree-Fock-Bogoliubov
  *
  *    D is split into four quadrants
  *
  *    | alpha alpha | alpha beta |
  *    ---------------------------- 
  *    | beta alpha  | beta beta  | 
  *
  *    The ghf orbitals are a combination of alpha and beta spins. In this case
  *    the mixed spin blocks have a contribution too. The same spin blocks have 
  *    the form
  *
  *    D_{mu nu}^{alpha alpha} = sum_{lambda sigma}( mu lambda|| nu sigma)k_{lambda sigma}^{alpha alpha} 
  *
  *    The mixed spin block doesn't have a coulomb term but it has an exchange 
  *    term.  
  *
  *    D_{mu nu}^{alpha beta} = sum_{lambda sigma}( mu lambda| nu sigma)k_{lambda sigma}^{alpha beta} -
  *     ( mu sigma| nu lambda)k_{lambda sigma}^{beta alpha} 
  *
  *    Pass each quadrant to be contracted with the two electron integrals.
  *
  */

    Eigen::MatrixXd t_p ;

    D.setZero() ;

    /* Do the pairing terms for the alpha alpha block */
    pairing_field( intarr, k.block( 0, 0, nb, nb), D.block( 0, 0, nb, nb), nb) ;

    /* Do the twin terms for the alpha alpha block */
    twin_field( intarr, k.block( 0, 0, nb, nb), D.block( 0, 0, nb, nb), nb) ;

    /* Do the pairing terms for the beta beta block */
    pairing_field( intarr, k.block( nb, nb, nb, nb), D.block( nb, nb, nb, nb), nb) ;

    /* Do the twin terms for the beta beta block */
    twin_field( intarr, k.block( nb, nb, nb, nb), D.block( nb, nb, nb, nb), nb) ;

    /* Do the pairing terms for the alpha beta block */
    pairing_field( intarr, k.block( 0, nb, nb, nb), D.block( 0, nb, nb, nb), nb) ;

    /* Do the twin terms for the alpha beta block */
    twin_field( intarr, k.block( nb, 0, nb, nb), D.block( 0, nb, nb, nb), nb) ;

    /* Do the pairing terms for the beta alpha block */
    pairing_field( intarr, k.block( nb, 0, nb, nb), D.block( nb, 0, nb, nb), nb) ;

    /* Do the twin terms for the beta alpha block */
    twin_field( intarr, k.block( 0, nb, nb, nb), D.block( nb, 0, nb, nb), nb) ;

    return ;
} ;

  void ctrPairg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> k, 
    Eigen::Ref<Eigen::MatrixXcd> D, const int nb) {

 /* 
  *    This routine contracts the two electron integrals to build the pairing
  *    field for Hartree-Fock-Bogoliubov
  *
  *    D is split into four quadrants
  *
  *    | alpha alpha | alpha beta |
  *    ---------------------------- 
  *    | beta alpha  | beta beta  | 
  *
  *    The ghf orbitals are a combination of alpha and beta spins. In this case
  *    the mixed spin blocks have a contribution too. The same spin blocks have 
  *    the form
  *
  *    D_{mu nu}^{alpha alpha} = sum_{lambda sigma}( mu lambda|| nu sigma)k_{lambda sigma}^{alpha alpha} 
  *
  *    The mixed spin block doesn't have a coulomb term but it has an exchange 
  *    term.  
  *
  *    D_{mu nu}^{alpha beta} = sum_{lambda sigma}( mu lambda| nu sigma)k_{lambda sigma}^{alpha beta} -
  *     ( mu sigma| nu lambda)k_{lambda sigma}^{beta alpha} 
  *
  *    Pass each quadrant to be contracted with the two electron integrals.
  *
  */

    D.setZero() ;

    /* Do the pairing terms for the alpha alpha block */
    pairing_field( intarr, k.block( 0, 0, nb, nb), D.block( 0, 0, nb, nb), nb) ;

    /* Do the twin terms for the alpha alpha block */
    twin_field( intarr, k.block( 0, 0, nb, nb), D.block( 0, 0, nb, nb), nb) ;

    /* Do the pairing terms for the beta beta block */
    pairing_field( intarr, k.block( nb, nb, nb, nb), D.block( nb, nb, nb, nb), nb) ;

    /* Do the twin terms for the beta beta block */
    twin_field( intarr, k.block( nb, nb, nb, nb), D.block( nb, nb, nb, nb), nb) ;

    /* Do the pairing terms for the alpha beta block */
    pairing_field( intarr, k.block( 0, nb, nb, nb), D.block( 0, nb, nb, nb), nb) ;

    /* Do the twin terms for the alpha beta block */
    twin_field( intarr, k.block( nb, 0, nb, nb), D.block( 0, nb, nb, nb), nb) ;

    /* Do the pairing terms for the beta alpha block */
    pairing_field( intarr, k.block( nb, 0, nb, nb), D.block( nb, 0, nb, nb), nb) ;

    /* Do the twin terms for the beta alpha block */
    twin_field( intarr, k.block( 0, nb, nb, nb), D.block( nb, 0, nb, nb), nb) ;

    return ;

} ;

void ctrPairr( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcd> k, Eigen::Ref<Eigen::MatrixXcd> D, const int nb) {

/* 
  This routine contracts the two electron integrals to build the pairing
  field for Hartree-Fock-Bogoliubov for the restricted model.
*/

    Eigen::MatrixXcd t ( nb, nb) ;
    t = -k ;
    D.setZero() ;

    /* Do the pairing terms for the alpha beta block */
    pairing_field( intarr, k, D, nb) ;

    /* Do the twin terms for the alpha beta block */
    twin_field( intarr, t, D, nb) ;

    t.resize( 0, 0) ;

    return ;

} ;

void ctrPairr( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXd> k, Eigen::Ref<Eigen::MatrixXd> D, const int nb) {

/* 
  This routine contracts the two electron integrals to build the pairing
  field for Hartree-Fock-Bogoliubov for the restricted model.
*/
    Eigen::MatrixXd t ( nb, nb) ;
    t = -k ;
    D.setZero() ;

    /* Do the pairing terms for the alpha beta block */
    pairing_field( intarr, k, D, nb) ;

    /* Do the twin terms for the alpha beta block */
    twin_field( intarr, t, D, nb) ;

    return ;

} ;

  int coulblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> G, const int nbasis) {
    /* Given a spin density block contract the coulomb terms into G 
     Contract coulomb integrals for complex matrices */
    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for ( int t = 0; t < n2ei; t++) {
      /*  ( 1 1| 2 2)
       *  ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;

      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i , i) += p( i, i)*val ;

      } else if ( ! ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i , j) += p( j, j)*val ;
        G( j , i) += p( j, j)*val ;
        G( j , j) += (p( j, i) + p( i, j))*val ;

      } else if ( ! ieqj && ! keql && i == k && j == l ){

        /* ( i j| i j) */

        G( i, j) += (p( i, j) + p( j, i))*val ;
        G( j, i) += (p( i, j) + p( j, i))*val ;

      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, i) += p( k, k)*val ;
        G( k, k) += p( i, i)*val ;

      } else if ( ieqj && ! keql && i == k ){

        /* ( i i| i l) */

        G( l, i) += p( i, i)*val ;
        G( i, l) += p( i, i)*val ;
        G( i, i) += (p( i, l) + p( l, i))*val ;

      } else if ( ! ieqj && ! keql && j == l && i != k){

        /* ( i j| k j) */

        G( i, j) += (p( k, j) + p( j, k))*val ;
        G( j, i) += (p( k, j) + p( j, k))*val ;
        G( k, j) += (p( i, j) + p( j, i))*val ;
        G( j, k) += (p( i, j) + p( j, i))*val ;

      } else if ( ! ieqj && keql && i != k && j != k) {

        /* ( i j| k k) */

        G( i, j) += p( k, k)*val ;
        G( j, i) += p( k, k)*val ;
        G( k, k) += (p( j, i) + p( i, j))*val ;

      } else if ( ! ieqj && ! keql && j == k && i != l ) {

        /* ( i j| j l) */

        G( i, j)+= (p( j, l) + p( l, j))*val ;
        G( j, i)+= (p( j, l) + p( l, j))*val ;
        G( j, l)+= (p( i, j) + p( j, i))*val ;
        G( l, j)+= (p( i, j) + p( j, i))*val ;

      } else if ( ! ieqj && ! keql && i == k && j != l) {

        /* ( i j| i l) */

        G( i, j)+= (p( i, l) + p( l, i))*val ;
        G( j, i)+= (p( i, l) + p( l, i))*val ;
        G( i, l)+= (p( i, j) + p( j, i))*val ;
        G( l, i)+= (p( i, j) + p( j, i))*val ;

      } else if ( ieqj && ! keql && i != k && i != l) {

        /* ( i i| k l) */

        G( i, i) += (p( k, l) + p( l, k))*val ;
        G( k, l) += p( i, i)*val ;
        G( l, k) += p( i, i)*val ;

      } else {

        /* ( i j| k l) */

        G( i, j)+= (p( k, l) + p( l, k))*val ;
        G( j, i)+= (p( k, l) + p( l, k))*val ;
        G( k, l)+= (p( i, j) + p( j, i))*val ;
        G( l, k)+= (p( i, j) + p( j, i))*val ;

      }


      }


    return 0 ;

  } ;

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> G, const int nbasis) {

    /* Given a spin density block contract the exchange terms into G 
     Contract exchange integrals for real matrices */

    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
      
      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i, i) += -p( i, i)*val ;

      } else if ( ! ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i, j) += -p( j, j)*val ;
        G( j, i) += -p( j, j)*val ;
        G( j, j) += -(p( j, i) + p( i, j))*val ;

      } else if ( ! ieqj && ! keql && i == k && j == l){

        /* ( i j| i j) */

        G( i, j) += -p( j, i)*val ;
        G( j, i) += -p( i, j)*val ;
        G( j, j) += -p( i, i)*val ;
        G( i, i) += -p( j, j)*val ;


      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, k) += -p( i, k)*val ;
        G( k, i) += -p( k, i)*val ;

      } else if ( ieqj && ! keql && i == k) {

        /* ( i i| i l) */

        G( i, i) += -(p( i, l) + p( l, i))*val ;
        G( l, i) += -p( i, i)*val ;
        G( i, l) += -p( i, i)*val ;

      } else if ( ! ieqj && ! keql && i != k && j == l ) {

        /* ( i j| k j) */

        G( j, j) += -(p( i, k) + p( k, i))*val ;
        G( i, j) += -p( j, k)*val ;
        G( j, k) += -p( i, j)*val ;
        G( i, k) += -p( j, j)*val ;
        G( k, j) += -p( j, i)*val ;
        G( j, i) += -p( k, j)*val ;
        G( k, i) += -p( j, j)*val ;

      } else if ( ! ieqj && keql && i != k && j != l ) {

        /* ( i j| k k) */

        G( i, k) += -p( j, k)*val ;
        G( j, k) += -p( i, k)*val ;
        G( k, j) += -p( k, i)*val ;
        G( k, i) += -p( k, j)*val ;

      } else if ( ! ieqj && ! keql && i != l && j == k ) {

        /* ( i j| j l) */

        G( j, j) += -(p( l, i) + p( i, l))*val ;
        G( i, l) += -p( j, j)*val ;
        G( i, j) += -p( j, l)*val ;
        G( j, l) += -p( i, j)*val ;
        G( l, j) += -p( j, i)*val ;
        G( j, i) += -p( l, j)*val ;
        G( l, i) += -p( j, j)*val ;

      } else if ( ! ieqj && ! keql && i == k && j != l ) {

        /* ( i j| i l) */

        G( i, i) += -(p( l, j) + p( j, l))*val ;
        G( i, l) += -p( j, i)*val ;
        G( j, l) += -p( i, i)*val ;
        G( j, i) += -p( i, l)*val ;
        G( i, j) += -p( l, i)*val ;
        G( l, j) += -p( i, i)*val ;
        G( l, i) += -p( i, j)*val ;

      } else if ( ieqj && ! keql && i != k && i != l ) {

        /* ( i i| k l) */

        G( i, l) += -p( i, k)*val ;
        G( i, k) += -p( i, l)*val ;
        G( k, i) += -p( l, i)*val ;
        G( l, i) += -p( k, i)*val ;

      } else {

        /* ( i j| k l) */

        G( i, l) += -p( j, k)*val ;
        G( i, k) += -p( j, l)*val ;
        G( j, l) += -p( i, k)*val ;
        G( j, k) += -p( i, l)*val ;
        G( k, j) += -p( l, i)*val ;
        G( l, j) += -p( k, i)*val ;
        G( k, i) += -p( l, j)*val ;
        G( l, i) += -p( k, j)*val ;

      }

      }

    return 0 ;

  } ;

  int coulblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> G, const int nbasis) {

    /* Given a spin density block contract the coulomb terms into G 
     Contract coulomb integrals for real matrices */

    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
      
      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i , i) += p( i, i)*val ;

      } else if ( ! ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i , j) += p( j, j)*val ;
        G( j , i) += p( j, j)*val ;
        G( j , j) += (p( j, i) + p( i, j))*val ;

      } else if ( ! ieqj && i == k && j == l ){

        /* ( i j| i J) */

        G( i, j) += (p( j, i) + p( i, j))*val ;
        G( j, i) += (p( i, j) + p( j, i))*val ;

      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, i) += p( k, k)*val ;
        G( k, k) += p( i, i)*val ;

      } else if ( ieqj && ! keql && i == k ){

        /* ( i i| i l) */

        G( i, i) += (p( i, l) + p( l, i))*val ;
        G( l, i) += p( i, i)*val ;
        G( i, l) += p( i, i)*val ;

      } else if ( ! ieqj && ! keql && j == l && i != k){

        /* ( i j| k j) */

        G( i, j) += (p( k, j) + p( j, k))*val ;
        G( j, i) += (p( k, j) + p( j, k))*val ;
        G( k, j) += (p( i, j) + p( j, i))*val ;
        G( j, k) += (p( i, j) + p( j, i))*val ;

      } else if ( ! ieqj && keql && i != k && j != k) {

        /* ( i j| k k) */

        G( i, j) += p( k, k)*val ;
        G( j, i) += p( k, k)*val ;
        G( k, k) += (p( j, i) + p( i, j))*val ;

      } else if ( ieqj && ! keql && i != k && i != l) {

        /* ( i i| k l) */

        G( i, i) += (p( k, l) + p( l, k))*val ;
        G( k, l) += p( i, i)*val ;
        G( l, k) += p( i, i)*val ;

      } else if ( ! ieqj && ! keql && j == k ) {

        /* ( i j| j l) */

        G( i, j)+= (p( j, l) + p( l, j))*val ;
        G( j, i)+= (p( j, l) + p( l, j))*val ;
        G( j, l)+= (p( i, j) + p( j, i))*val ;
        G( l, j)+= (p( i, j) + p( j, i))*val ;

      } else if ( ! ieqj && ! keql && i == k && j != l) {

        /* ( i j| i l) */

        G( i, j)+= (p( i, l) + p( l, i))*val ;
        G( j, i)+= (p( i, l) + p( l, i))*val ;
        G( i, l)+= (p( i, j) + p( j, i))*val ;
        G( l, i)+= (p( i, j) + p( j, i))*val ;

      } else {

        /* ( i j| k l) */

        G( i, j)+= (p( k, l) + p( l, k))*val ;
        G( j, i)+= (p( k, l) + p( l, k))*val ;
        G( k, l)+= (p( i, j) + p( j, i))*val ;
        G( l, k)+= (p( i, j) + p( j, i))*val ;

      }

      }

    return 0 ;

  } ;

//  This is an implementation of the exchange block contraction which I know works
//  I implemented a new version to be consistent with the other routines and I want
//  to preserve this until I'm confident my new routine has zero bugs.  
//  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
//    Eigen::Ref<Eigen::MatrixXd> G, const int nbasis) {
//
//    /* Given a spin density block contract the exchange terms into G 
//     Contract exchange integrals for real matrices */
//
//    int ntt ;
//    int n2ei ;
//    int i ;
//    int j ;
//    int k ;
//    int l ;
//    double val ;
//    bool ieqj ;
//    bool keql ;
//
//    ntt = nbasis*(nbasis + 1)/2 ;
//    n2ei = ntt*(ntt+1)/2 ;
//
//    for(int t = 0; t < n2ei; t++ ) {
//      /*  ( 1 1| 2 2)
// *        ( i j| k l) = val */
//      i = intarr[t].r_i() ;
//      j = intarr[t].r_j() ;
//      k = intarr[t].r_k() ;
//      l = intarr[t].r_l() ;
//      val = intarr[t].r_v() ;
////   
//      ieqj = ( ( i == j) ? true : false) ;
//      keql = ( ( k == l) ? true : false) ;
//      
//      if ( ieqj && keql && i == k ){
//
//        /* ( i i| i i) */
//
//        G( i , i) += -p( i, i)*val ;
//
//      } else if ( not ieqj && keql && j == k){
//
//        /* ( i j| j j) */
//
//        G( i , j) += -p( j, j)*val ;
//        G( j , i) += -p( j, j)*val ;
//        G( j , j) += -(p( j, i) + p( i, j))*val ;
//
//      } else if ( not ieqj && i == k && j == l){
//
//        /* ( i j| i j) */
//
//        G( i, j) += -p( i, j)*val ;
//        G( j, i) += -p( j, i)*val ;
//        G( i, i) += -p( j, j)*val ;
//        G( j, j) += -p( i, i)*val ;
//
//      } else if ( ieqj && keql && i != k ){
//
//        /* ( i i| k k) */
//
//        G( i, k) += -p( k, i)*val ;
//        G( k, i) += -p( i, k)*val ;
//
//      } else if ( ieqj && not keql && i == k ){
//
//        /* ( i i| i l) */
//
//        G( i, i) += -(p( i, l) + p( l, i))*val ;
//        G( l, i) += -p( i, i)*val ;
//        G( i, l) += -p( i, i)*val ;
//
//      } else if ( not ieqj && not keql && j == l && i != k){
//
//        /* ( i j| k j) */
//
//        G( i, j) += -p( k, j)*val ;
//        G( j, i) += -p( j, k)*val ;
//        G( k, j) += -p( i, j)*val ;
//        G( j, k) += -p( j, i)*val ;
//        G( i, k) += -p( j, j)*val ;
//        G( k, i) += -p( j, j)*val ;
//        G( j, j) += -(p( i, k) + p( k, i))*val ;
//
//      } else if ( not ieqj && keql && i != k ) {
//
//        /* ( i j| k k) */
//
//        G( i, k) += -p( k, j)*val ;
//        G( j, k) += -p( k, i)*val ;
//        G( k, j) += -p( i, k)*val ;
//        G( k, i) += -p( j, k)*val ;
//
//      } else if ( ieqj && not keql && i != k && i != l) {
//
//        /* ( i i| k l) */
//
//        G( i, l) += -p( k, i)*val ;
//        G( i, k) += -p( l, i)*val ;
//        G( k, i) += -p( i, l)*val ;
//        G( l, i) += -p( i, k)*val ;
//
//      } else if ( not ieqj && not keql && j == k ) {
//
//        /* ( i j| j l) */
//
//        G( i, l) += -p( j, j)*val ;
//        G( l, i) += -p( j, j)*val ;
//        G( l, j) += -p( i, j)*val ;
//        G( j, l) += -p( j, i)*val ;
//        G( i, j) += -p( l, j)*val ;
//        G( j, i) += -p( j, l)*val ;
//        G( j, j) += -(p( i, l) + p( l, i))*val ;
//
//      } else if ( not ieqj && not keql && i == k && j != l) {
//
//        /* ( i j| i l) */
//
//        G( i, l) += -p( i, j)*val ;
//        G( l, i) += -p( j, i)*val ;
//        G( l, j) += -p( i, i)*val ;
//        G( j, l) += -p( i, i)*val ;
//        G( i, j) += -p( i, l)*val ;
//        G( j, i) += -p( l, i)*val ;
//        G( i, i) += -(p( j, l) + p( l, j))*val ;
//
//      } else {
//
//        /* ( i j| k l) */
//
//        G( i, l)+= -p( k, j)*val ;
//        G( j, l)+= -p( k, i)*val ;
//        G( k, j)+= -p( i, l)*val ;
//        G( j, k)+= -p( l, i)*val ;
//        G( i, k)+= -p( l, j)*val ;
//        G( l, j)+= -p( i, k)*val ;
//        G( k, i)+= -p( j, l)*val ;
//        G( l, i)+= -p( j, k)*val ;
//
//      }
//
//      }
//
//    return 0 ;
//
//  } ; 

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> G, const int nbasis) {

    /* Given a spin density block contract the exchange terms into G 
     Contract exchange integrals for real matrices */

    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
      
      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i, i) += -p( i, i)*val ;

      } else if ( not ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i, j) += -p( j, j)*val ;
        G( j, i) += -p( j, j)*val ;
        G( j, j) += -(p( j, i) + p( i, j))*val ;

      } else if ( not ieqj && not keql && i == k && j == l){

        /* ( i j| i j) */

        G( i, j) += -p( j, i)*val ;
        G( j, i) += -p( i, j)*val ;
        G( j, j) += -p( i, i)*val ;
        G( i, i) += -p( j, j)*val ;


      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, k) += -p( i, k)*val ;
        G( k, i) += -p( k, i)*val ;

      } else if ( ieqj && not keql && i == k) {

        /* ( i i| i l) */

        G( i, i) += -(p( i, l) + p( l, i))*val ;
        G( l, i) += -p( i, i)*val ;
        G( i, l) += -p( i, i)*val ;

      } else if ( not ieqj && not keql && i != k && j == l ) {

        /* ( i j| k j) */

        G( j, j) += -(p( i, k) + p( k, i))*val ;
        G( i, j) += -p( j, k)*val ;
        G( j, k) += -p( i, j)*val ;
        G( i, k) += -p( j, j)*val ;
        G( k, j) += -p( j, i)*val ;
        G( j, i) += -p( k, j)*val ;
        G( k, i) += -p( j, j)*val ;

      } else if ( not ieqj && keql && i != k && j != l ) {

        /* ( i j| k k) */

        G( i, k) += -p( j, k)*val ;
        G( j, k) += -p( i, k)*val ;
        G( k, j) += -p( k, i)*val ;
        G( k, i) += -p( k, j)*val ;

      } else if ( not ieqj && not keql && i != l && j == k ) {

        /* ( i j| j l) */

        G( j, j) += -(p( l, i) + p( i, l))*val ;
        G( i, l) += -p( j, j)*val ;
        G( i, j) += -p( j, l)*val ;
        G( j, l) += -p( i, j)*val ;
        G( l, j) += -p( j, i)*val ;
        G( j, i) += -p( l, j)*val ;
        G( l, i) += -p( j, j)*val ;

      } else if ( not ieqj && not keql && i == k && j != l ) {

        /* ( i j| i l) */

        G( i, i) += -(p( l, j) + p( j, l))*val ;
        G( i, l) += -p( j, i)*val ;
        G( j, l) += -p( i, i)*val ;
        G( j, i) += -p( i, l)*val ;
        G( i, j) += -p( l, i)*val ;
        G( l, j) += -p( i, i)*val ;
        G( l, i) += -p( i, j)*val ;

      } else if ( ieqj && not keql && i != k && i != l ) {

        /* ( i i| k l) */

        G( i, l) += -p( i, k)*val ;
        G( i, k) += -p( i, l)*val ;
        G( k, i) += -p( l, i)*val ;
        G( l, i) += -p( k, i)*val ;

      } else {

        /* ( i j| k l) */

        G( i, l) += -p( j, k)*val ;
        G( i, k) += -p( j, l)*val ;
        G( j, l) += -p( i, k)*val ;
        G( j, k) += -p( i, l)*val ;
        G( k, j) += -p( l, i)*val ;
        G( l, j) += -p( k, i)*val ;
        G( k, i) += -p( l, j)*val ;
        G( l, i) += -p( k, j)*val ;

      }

      }

    return 0 ;

  } ;

  int pairing_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> D, const int nbasis) {
    /* Evaluate the pairing field for HFB type wavefunctions */
    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
//
      if ( ieqj && keql && i == k) {

        /* ( i i| i i)*/

        D(i,i) += val*p(i,i) ;

      } else if ( not ieqj && keql && j == k) {

        /* ( i j| j j)*/

        D(i,j) += val*p(j,j) ;
        D(j,i) += val*p(j,j) ;
        D(j,j) += val*(p(i,j) + p(j,i)) ;

      } else if ( not ieqj && not keql && i == k && j ==l) {

        /* ( i j| i j)*/

        D(i,i) += val*p(j,j) ;
        D(j,j) += val*p(i,i) ;
        D(i,j) += val*p(j,i) ;
        D(j,i) += val*p(i,j) ;

      } else if ( ieqj && keql && i !=k) {

        /* ( i i| k k) */

        D(i,k) += val*p(i,k) ;
        D(k,i) += val*p(k,i) ;

      } else if ( ieqj && not keql && i == k){

        /* ( i i| i l) */

        D(i,l) += val*p(i,i) ;
        D(l,i) += val*p(i,i) ;
        D(i,i) += val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i != k && j == l){

        /* ( i j| k j) */

        D(i,k) += val*p(j,j) ;
        D(j,k) += val*p(i,j) ;
        D(i,j) += val*p(j,k) ;
        D(k,i) += val*p(j,j) ;
        D(j,i) += val*p(k,j) ;
        D(k,j) += val*p(j,i) ;
        D(j,j) += val*(p(i,k) + p(k,i)) ;

      } else if ( not ieqj && keql && i != k && j != l) {

        /* ( i j| k k) */

        D(i,k) += val*p(j,k) ;
        D(k,i) += val*p(k,j) ;
        D(j,k) += val*p(i,k) ;
        D(k,j) += val*p(k,i) ;

      } else if ( not ieqj && not keql && i!=l && j == k) {

        /* ( i j| j l) */

        D(i,j) += val*p(j,l) ;
        D(i,l) += val*p(j,j) ;
        D(j,l) += val*p(i,j) ;
        D(j,i) += val*p(l,j) ;
        D(l,i) += val*p(j,j) ;
        D(l,j) += val*p(j,i) ;
        D(j,j) += val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i == k && j != l){

        /* ( i j| i l) */

        D(i,l) += val*p(j,i) ;
        D(j,i) += val*p(i,l) ;
        D(j,l) += val*p(i,i) ;
        D(l,j) += val*p(i,i) ;
        D(i,j) += val*p(l,i) ;
        D(l,i) += val*p(i,j) ;
        D(i,i) += val*(p(j,l) + p(l,j)) ;

      } else if ( ieqj && not keql && i != k && i != l){

        /* ( i i| k l) */

        D(i,k) += val*p(i,l) ;
        D(i,l) += val*p(i,k) ;
        D(k,i) += val*p(l,i) ;
        D(l,i) += val*p(k,i) ;

      } else {

        /* ( i j| k l) */

        D(i,k) += val*p(j,l) ;
        D(i,l) += val*p(j,k) ;
        D(j,k) += val*p(i,l) ;
        D(j,l) += val*p(i,k) ;
        D(k,i) += val*p(l,j) ;
        D(l,i) += val*p(k,j) ;
        D(l,j) += val*p(k,i) ;
        D(k,j) += val*p(l,i) ;

      }

      }

  return 0 ;

  } ;

  int twin_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> D, const int nbasis) {
    /* Evaluate the twin field for HFB type wavefunctions */
    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
//
      if ( ieqj && keql && i == k) {

        /* ( i i| i i)*/

        D(i,i) += -val*p(i,i) ;

      } else if ( not ieqj && keql && j == k) {

        /* ( i j| j j)*/

        D(i,j) += -val*p(j,j) ;
        D(j,i) += -val*p(j,j) ;
        D(j,j) += -val*(p(i,j) + p(j,i)) ;

      } else if ( ! ieqj && ! keql && i == k && j == l) {

        /* ( i j| i j)*/

        D(i,i) += -val*p(j,j) ;
        D(j,j) += -val*p(i,i) ;
        D(i,j) += -val*p(i,j) ;
        D(j,i) += -val*p(j,i) ;

      } else if ( ieqj && keql && i != k) {

        /* ( i i| k k) */

        D(i,k) += -val*p(k,i) ;
        D(k,i) += -val*p(i,k) ;

      } else if ( ieqj && not keql && i == k){

        /* ( i i| i l) */

        D(i,l) += -val*p(i,i) ;
        D(l,i) += -val*p(i,i) ;
        D(i,i) += -val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i != k && j == l){

        /* ( i j| k j) */

        D(i,k) += -val*p(j,j) ;
        D(k,i) += -val*p(j,j) ;
        D(j,k) += -val*p(j,i) ;
        D(i,j) += -val*p(k,j) ;
        D(j,i) += -val*p(j,k) ;
        D(k,j) += -val*p(i,j) ;
        D(j,j) += -val*(p(i,k) + p(k,i)) ;

      } else if ( not ieqj && keql && i != k && j != l) {

        /* ( i j| k k) */

        D(i,k) += -val*p(k,j) ;
        D(j,k) += -val*p(k,i) ;
        D(k,i) += -val*p(j,k) ;
        D(k,j) += -val*p(i,k) ;

      } else if ( not ieqj && not keql && i!=l && j == k) {

        /* ( i j| j l) */

        D( i, j) += -val*p( l, j) ;
        D( j, l) += -val*p( j, i) ;
        D( i, l) += -val*p( j, j) ;
        D( j, i) += -val*p( j, l) ;
        D( l, i) += -val*p( j, j) ;
        D( l, j) += -val*p( i, j) ;
        D( j, j) += -val*(p( i, l) + p( l, i)) ;

      } else if ( not ieqj && not keql && i == k && j != l){

        /* ( i j| i l) */

        D( j, i) += -val*p( l, i) ;
        D( j, l) += -val*p( i, i) ;
        D( i, l) += -val*p( i, j) ;
        D( i, j) += -val*p( i, l) ;
        D( l, i) += -val*p( j, i) ;
        D( l, j) += -val*p( i, i) ;
        D( i, i) += -val*(p( j, l) + p( l, j)) ;

      } else if ( ieqj && not keql && i != k && i != l){

        /* ( i i| k l) */

        D( i, k) += -val*p( l, i) ;
        D( i, l) += -val*p( i, k) ;
        D( k, i) += -val*p( i, l) ;
        D( l, i) += -val*p( i, k) ;

      } else {

        /* ( i j| k l) */

        D( i, k) += -val*p( l, j) ;
        D( j, k) += -val*p( l, i) ;
        D( j, l) += -val*p( k, i) ;
        D( i, l) += -val*p( j, k) ;
        D( k, i) += -val*p( j, l) ;
        D( l, i) += -val*p( j, k) ;
        D( l, j) += -val*p( i, k) ;
        D( k, j) += -val*p( i, l) ;

      }

      }

  return 0 ;

  } ;

  int pairing_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcd> p, 
    Eigen::Ref<Eigen::MatrixXcd> D, const int nbasis) {
    /* Evaluate the pairing field for HFB type wavefunctions */
    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
//
      if ( ieqj && keql && i == k) {

        /* ( i i| i i)*/

        D(i,i) += val*p(i,i) ;

      } else if ( not ieqj && keql && j == k) {

        /* ( i j| j j)*/

        D(i,j) += val*p(j,j) ;
        D(j,i) += val*p(j,j) ;
        D(j,j) += val*(p(i,j) + p(j,i)) ;

      } else if ( not ieqj && not keql && i == k && j ==l) {

        /* ( i j| i j)*/

        D(i,i) += val*p(j,j) ;
        D(j,j) += val*p(i,i) ;
        D(i,j) += val*p(j,i) ;
        D(j,i) += val*p(i,j) ;

      } else if ( ieqj && keql && i !=k) {

        /* ( i i| k k) */

        D(i,k) += val*p(i,k) ;
        D(k,i) += val*p(k,i) ;

      } else if ( ieqj && not keql && i == k){

        /* ( i i| i l) */

        D(i,l) += val*p(i,i) ;
        D(l,i) += val*p(i,i) ;
        D(i,i) += val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i != k && j == l){

        /* ( i j| k j) */

        D(i,k) += val*p(j,j) ;
        D(j,k) += val*p(i,j) ;
        D(i,j) += val*p(j,k) ;
        D(k,i) += val*p(j,j) ;
        D(j,i) += val*p(k,j) ;
        D(k,j) += val*p(j,i) ;
        D(j,j) += val*(p(i,k) + p(k,i)) ;

      } else if ( not ieqj && keql && i != k && j != l) {

        /* ( i j| k k) */

        D(i,k) += val*p(j,k) ;
        D(k,i) += val*p(k,j) ;
        D(j,k) += val*p(i,k) ;
        D(k,j) += val*p(k,i) ;

      } else if ( not ieqj && not keql && i!=l && j == k) {

        /* ( i j| j l) */

        D(i,j) += val*p(j,l) ;
        D(i,l) += val*p(j,j) ;
        D(j,l) += val*p(i,j) ;
        D(j,i) += val*p(l,j) ;
        D(l,i) += val*p(j,j) ;
        D(l,j) += val*p(j,i) ;
        D(j,j) += val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i == k && j != l){

        /* ( i j| i l) */

        D(i,l) += val*p(j,i) ;
        D(j,i) += val*p(i,l) ;
        D(j,l) += val*p(i,i) ;
        D(l,j) += val*p(i,i) ;
        D(i,j) += val*p(l,i) ;
        D(l,i) += val*p(i,j) ;
        D(i,i) += val*(p(j,l) + p(l,j)) ;

      } else if ( ieqj && not keql && i != k && i != l){

        /* ( i i| k l) */

        D(i,k) += val*p(i,l) ;
        D(i,l) += val*p(i,k) ;
        D(k,i) += val*p(l,i) ;
        D(l,i) += val*p(k,i) ;

      } else {

        /* ( i j| k l) */

        D(i,k) += val*p(j,l) ;
        D(i,l) += val*p(j,k) ;
        D(j,k) += val*p(i,l) ;
        D(j,l) += val*p(i,k) ;
        D(k,i) += val*p(l,j) ;
        D(l,i) += val*p(k,j) ;
        D(l,j) += val*p(k,i) ;
        D(k,j) += val*p(l,i) ;

      }

      }

  return 0 ;

  } ;

  int twin_field( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXd> p, 
    Eigen::Ref<Eigen::MatrixXd> D, const int nbasis) {
    /* Evaluate the twin field for HFB type wavefunctions */
    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    double val ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*
          ( 1 1| 2 2)
          ( i j| k l) = val
      */
      i = intarr[t].r_i() ;
      j = intarr[t].r_j() ;
      k = intarr[t].r_k() ;
      l = intarr[t].r_l() ;
      val = intarr[t].r_v() ;
//
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
//
      if ( ieqj && keql && i == k) {

        /* ( i i| i i)*/

        D(i,i) += -val*p(i,i) ;

      } else if ( not ieqj && keql && j == k) {

        /* ( i j| j j)*/

        D(i,j) += -val*p(j,j) ;
        D(j,i) += -val*p(j,j) ;
        D(j,j) += -val*(p(i,j) + p(j,i)) ;

      } else if ( not ieqj && not keql && i == k && j ==l) {

        /* ( i j| i j)*/

        D(i,i) += -val*p(j,j) ;
        D(j,j) += -val*p(i,i) ;
        D(i,j) += -val*p(i,j) ;
        D(j,i) += -val*p(j,i) ;

      } else if ( ieqj && keql && i !=k) {

        /* ( i i| k k) */

        D(i,k) += -val*p(k,i) ;
        D(k,i) += -val*p(i,k) ;

      } else if ( ieqj && not keql && i == k){

        /* ( i i| i l) */

        D(i,l) += -val*p(i,i) ;
        D(l,i) += -val*p(i,i) ;
        D(i,i) += -val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i != k && j == l){

        /* ( i j| k j) */

        D(i,k) += -val*p(j,j) ;
        D(k,i) += -val*p(j,j) ;
        D(j,k) += -val*p(j,i) ;
        D(i,j) += -val*p(k,j) ;
        D(j,i) += -val*p(j,k) ;
        D(k,j) += -val*p(i,j) ;
        D(j,j) += -val*(p(i,k) + p(k,i)) ;

      } else if ( not ieqj && keql && i != k && j != l) {

        /* ( i j| k k) */

        D(i,k) += -val*p(k,j) ;
        D(j,k) += -val*p(k,i) ;
        D(k,i) += -val*p(j,k) ;
        D(k,j) += -val*p(i,k) ;

      } else if ( not ieqj && not keql && i!=l && j == k) {

        /* ( i j| j l) */

        D(i,j) += -val*p(l,j) ;
        D(j,l) += -val*p(j,i) ;
        D(i,l) += -val*p(j,j) ;
        D(j,i) += -val*p(j,l) ;
        D(l,i) += -val*p(j,j) ;
        D(l,j) += -val*p(i,j) ;
        D(j,j) += -val*(p(i,l) + p(l,i)) ;

      } else if ( not ieqj && not keql && i == k && j != l){

        /* ( i j| i l) */

        D(j,i) += -val*p(l,i) ;
        D(j,l) += -val*p(i,i) ;
        D(i,l) += -val*p(i,j) ;
        D(i,j) += -val*p(i,l) ;
        D(l,i) += -val*p(j,i) ;
        D(l,j) += -val*p(i,i) ;
        D(i,i) += -val*(p(j,l) + p(l,j)) ;

      } else if ( ieqj && not keql && i != k && i != l){

        /* ( i i| k l) */

        D(i,k) += -val*p(l,i) ;
        D(i,l) += -val*p(i,k) ;
        D(k,i) += -val*p(i,l) ;
        D(l,i) += -val*p(i,k) ;

      } else {

        /* ( i j| k l) */

        D(i,k) += -val*p(l,j) ;
        D(j,k) += -val*p(l,i) ;
        D(j,l) += -val*p(k,i) ;
        D(i,l) += -val*p(j,k) ;
        D(k,i) += -val*p(j,l) ;
        D(l,i) += -val*p(j,k) ;
        D(l,j) += -val*p(i,k) ;
        D(k,j) += -val*p(i,l) ;

      }

      }

  return 0 ;

  } ;

/* Evaluate elements between Slater Determinants */
double tranden1 ( int& nele, int& nbasis, Eigen::Ref<Eigen::MatrixXd> wfna, Eigen::Ref<Eigen::MatrixXd> wfnb, Eigen::Ref<Eigen::MatrixXd> dabmat, Eigen::Ref<Eigen::MatrixXd> scr) {

/*
 * Given mo coeffcients for Slater Determinants, return the transition density
 * and the overlap for a single block.
 *
 * The determinants need to be in oao basis.
 *
 * sum_{kl} C_{vl}D(k | l)C_{ku}^{*}
 *
 * D =
 * [ d_{1_{a}1_{b}} d_{1_{a}2_{b}} d_{1_{a}3_{b}} ... ]
 * [ d_{2_{a}1_{b}} d_{2_{a}2_{b}} d_{2_{a}3_{b}} ... ]
 * [       |              |              |         |  ]
 * [ d_{n_{a}1_{b}} d_{n_{a}2_{b}} d_{n_{a}3_{b}} ... ]
 *
 *  d_{kl} = int psi_{k}^{*}(1)psi_{l}(1) 
 *
 * -This routine assumes both determinants have the same symmetry structure.
 * -This routine assumes that the determinants are in the orthongonal
 *    ao basis.
 * -This routine assumes that the determinants are normalized themselves
 *  so the transition density is implicity divided by 1.
 * -Returns a GHF style transition density matrix even if many elements
 *  are zero
 *
 * */

  double ovl ;

  /* Build D*/
  dabmat.block( 0, 0, nele, nele) = wfna.block( 0, 0, nbasis, nele).adjoint()*wfnb.block( 0, 0, nbasis, nele) ;

  /* Build D(k | l) */

  ovl = dabmat.block( 0, 0, nele, nele).determinant() ;

  /* Build the adjugate */
  scr = ovl*dabmat.block( 0, 0, nele, nele).inverse() ;

  dabmat = wfnb.block( 0, 0, nbasis, nele)*scr*wfna.block( 0, 0, nbasis, nele).transpose() ;

  return ovl ;

} ;

cd tranden1 ( int& nele, int& nbasis, Eigen::Ref<Eigen::MatrixXcd> wfna, Eigen::Ref<Eigen::MatrixXcd> wfnb, Eigen::Ref<Eigen::MatrixXcd> dabmat, Eigen::Ref<Eigen::MatrixXcd> scr) {

/*
 * Complex overloaded transition density routine
 * */

  cd ovl ;

  /* Build D*/
  dabmat.block( 0, 0, nele, nele) = wfna.block( 0, 0, nbasis, nele).adjoint()*wfnb.block( 0, 0, nbasis, nele) ;

  /* Build D(k | l) */

  ovl = dabmat.block( 0, 0, nele, nele).determinant() ;

  /* Build the adjugate */
  scr = ovl*dabmat.block( 0, 0, nele, nele).inverse() ;

  dabmat = wfnb.block( 0, 0, nbasis, nele)*scr*wfna.block( 0, 0, nbasis, nele).transpose() ;

  return ovl ;

} ;

double tranden2 ( int& ne1, int& ne2, int& nbasis, Eigen::Ref<Eigen::MatrixXd> wfna, Eigen::Ref<Eigen::MatrixXd> wfnb, Eigen::Ref<Eigen::MatrixXd> dabmat, Eigen::Ref<Eigen::MatrixXd> scr1) {

/*
 * Given mo coeffcients for an unrestricted determinant, find the transition density
 * for the alpha and beta components.
 *
 * The mos should be passed as alpha block/beta block.
 ^
 ^ I think this means we simply do tranden1 twice.
 * No sure if I can pass a block of an eigen matrix so i will use scratch for
 * now.
 *
 * */

  double ola, olb ;

  /* Build D*/
  ola = tranden1( ne1, nbasis, wfna.block( 0, 0, nbasis, nbasis), wfnb.block( 0, 0, nbasis, nbasis), dabmat.block( 0, 0, nbasis, nbasis), scr1.block( 0, 0, ne1, ne1)) ;

  olb = tranden1( ne2, nbasis, wfna.block( 0, nbasis, nbasis, nbasis), wfnb.block( 0, nbasis, nbasis, nbasis), dabmat.block( 0, nbasis, nbasis, nbasis), scr1.block( 0, 0, ne2, ne2)) ;

  return ola*olb ;

} ;

cd tranden2 ( int& ne1, int& ne2, int& nbasis, Eigen::Ref<Eigen::MatrixXcd> wfna, Eigen::Ref<Eigen::MatrixXcd> wfnb, Eigen::Ref<Eigen::MatrixXcd> dabmat, Eigen::Ref<Eigen::MatrixXcd> scr) {

  cd ola, olb ;

  /* Build D*/
  ola = tranden1( ne1, nbasis, wfna.block( 0, 0, nbasis, nbasis), wfnb.block( 0, 0, nbasis, nbasis), dabmat.block( 0, 0, nbasis, nbasis), scr.block( 0, 0, ne1, ne1)) ;

  olb = tranden1( ne2, nbasis, wfna.block( 0, nbasis, nbasis, nbasis), wfnb.block( 0, nbasis, nbasis, nbasis), dabmat.block( 0, nbasis, nbasis, nbasis), scr.block( 0, 0, ne2, ne2)) ;

  return ola*olb ;

} ;

//double obop ( common& com, Eigen::Ref<Eigen::MatrixXd> ouv, hfwfn& a, hfwfn& b) {
//
///*
// * One body operator.  Return <A|O|B> = sum_{ l k}<k|O|l> = sum{ u v} <u|O|v>C_{vl}D(k|l)C_{ku}^{*}
// * Passed two Slater determinants and matrix elements in an orthogonal basis, return the evaluated operator. */
//
//  double aob ;
//  Eigen::MatrixXd pvu ;
//  Eigen::MatrixXd omega ;
//
//  /* Build pvu */
// 
//  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.setZero() ;
//
// // tranden1( com, a, b, pvu) ;
//
//  omega = ouv*pvu ;
//  aob = omega.trace() ;
//
//  pvu.resize( 0, 0) ;
//  omega.resize( 0, 0) ;
//
//  return aob ;
//
//} ;
//
// cd obop ( common& com, Eigen::Ref<Eigen::MatrixXcd> ouv, hfwfn& a, hfwfn& b) {
//
///*
// * Overloaded for complex functions.
// */
//
//  cd aob ;
//  Eigen::MatrixXcd pvu ;
//  Eigen::MatrixXcd omega ;
//
//  /* Build pvu */
//
//  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.setZero() ;
//
//  //tranden ( com, a, b, pvu) ;
//
//  omega = ouv*pvu ;
//  aob = omega.trace() ;
//
//  return aob ;
//
//} ;
//
//template<typename >
//double fockop ( common& com, Eigen::Ref<Eigen::MatrixXd> h, std::vector<tei>& intarr, hfwfn& a,
//               hfwfn& b, double& ovl) { 
//
///*
// * Given a density matrix return 
// * */
//
//  double aob ;
//  Eigen::MatrixXd pvu ;
//  Eigen::MatrixXd omega ;
//  Eigen::MatrixXd f ;
//  Eigen::MatrixXd g ;
//
//  /* Build pvu */
//
//  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
//  f.resize( 2*com.nbas(), 2*com.nbas()) ;
//  g.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.setZero() ;
//  f.setZero() ;
//
////  ovl = tranden ( com, a, b, pvu) ;
//
//  ctr2eg( intarr, pvu, g, com.nbas()) ;
//
//  f =  h + g ;
//  g = h + f ;
//  omega = g*pvu ;
//  aob = 0.5*omega.trace() ;
//
//  g.resize( 0, 0) ;
//  f.resize( 0, 0) ;
//  pvu.resize( 0, 0) ;
//  omega.resize( 0, 0) ;
//
//  return aob ;
//
//} ;
//
// cd fockop ( common& com, Eigen::Ref<Eigen::MatrixXcd> h, std::vector<tei>& intarr, hfwfn& a, 
//                             hfwfn& b, cd& ovl) {
//
///*
//  While the fock operator is one body, it requires that we contract the two electron
//  integrals with the density to build it.  This wraps up that procedure.  This routine assumes
// 
//   - Everything has been put into an orthogonal ao basis.
//   - Regardless of the flavor of hf wfn, this treats everything as ghf meaning 2*nbas dimensions.
//*/
//
//  cd aob ;
//  cd pt5 = cd (0.5,0.0) ;
//  Eigen::MatrixXcd pvu ;
//  Eigen::MatrixXcd omega ;
//  Eigen::MatrixXcd f ;
//  Eigen::MatrixXcd g ;
//  Eigen::MatrixXcd mos ;
//
//  /* Build pvu */
//
//  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
//  f.resize( 2*com.nbas(), 2*com.nbas()) ;
//  g.resize( 2*com.nbas(), 2*com.nbas()) ;
//  mos.resize( 2*com.nbas(), 2*com.nbas()) ;
//  pvu.setZero() ;
//  f.setZero() ;
//  a.get_mos( mos) ;
//
////  ovl = tranden ( com, a, b, pvu) ;
//
//  ctr2eg( intarr, pvu, g, com.nbas()) ;
//
//  f =  h + g ;
//  g = h + f ;
//  omega = g*pvu ;
//  aob = pt5*omega.trace() ;
//
//  g.resize( 0, 0) ;
//  f.resize( 0, 0) ;
//  pvu.resize( 0, 0) ;
//  omega.resize( 0, 0) ;
//
//  return aob ;
//
//} ;
//
