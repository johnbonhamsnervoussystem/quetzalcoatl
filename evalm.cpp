#include "constants.h"
#include<iostream>
#include<vector>
#include<Eigen/Dense>
#include "tei.h"
#include "hfwfn.h"
#include "evalm.h"

/* evalm is a collection of routines that evaluates matrix elements  */

/* Given a vector of two electron integrals contract it with a density
 * and return the matrix. */

  int ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> g, const int nb) {

 /* G is split into four quadrants.  
  *
  *    | alpha alpha | alpha beta |
  *    ---------------------------- 
  *    | beta alpha  | beta beta  | 
  *
  *    For restricted Hartree-Fock calculations, the orbitals are either
  *    alpha or beta, not a combination of both so the mixed blocks are
  *    zero.  Further, the alpha and beta blocks are identical so we only
  *    need to calculate the coulomb and exchange for a single block.
  *
  *    G_{mu nu} = sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha + 
  *    p_{sigma lambda}^{beta beta })( mu nu| lambda sigma) - 
  *    sum_{lambda sigma}(p_{sigma lambda}^{alpha alpha)( mu nu| sigma lambda ) 
  *
  */

    Eigen::MatrixXf t_p ;

    g.setZero() ;
    t_p.resize( nb, nb) ;

    t_p = 2.0*p ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, t_p, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    t_p.resize( 0, 0) ;

    return 0 ;
} ;

  int ctr2er( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> g, const int nb) {

 /*  
  * This is an overloaded function with complex matrices for complex
  * restricted Hartree-Fock.
  */

    Eigen::MatrixXcf t_p ;

    g.setZero() ;
    t_p.resize( nb, nb) ;

    t_p = 2.0*p ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, t_p, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    t_p.resize( 0, 0) ;

    return 0 ;
} ;

  int ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> pt, 
     Eigen::Ref<Eigen::MatrixXf> p, Eigen::Ref<Eigen::MatrixXf> g, const int nb) {

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

    return 0 ;
} ;

  int ctr2eu( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> pt, 
     Eigen::Ref<Eigen::MatrixXcf> p, Eigen::Ref<Eigen::MatrixXcf> g, const int nb) {

 /*  
  * This is an overloaded function with complex matrices for complex
  * unrestricted Hartree-Fock.
  */
    g.setZero() ;

    /* Do the coulomb terms for the alpha alpha block  */
    coulblk( intarr, pt, g, nb) ;

    /* Do the exchange terms for the alpha alpha block  */
    exchblk( intarr, p, g, nb) ;

    return 0 ;
} ;

  int ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> g, const int nb) {

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

    Eigen::MatrixXf t_p ;

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

    return 0 ;
} ;

  int ctr2eg( std::vector<tei>& intarr, Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> g, const int nb) {

 /*   
  * This is the complex overloaded version of ctr2eg
  */
    Eigen::MatrixXcf t_p ;

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

    return 0 ;
} ;

  int coulblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> G, const int nbasis) {
    /* Given a spin density block contract the coulomb terms into G 
     Contract coulomb integrals for complex matrices */
    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    float val ;
    bool ieqk ;
    bool jeql ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
       *  ( i j| k l) = val */
      i = intarr[t].r_i() - 1 ;
      j = intarr[t].r_j() - 1 ;
      k = intarr[t].r_k() - 1 ;
      l = intarr[t].r_l() - 1 ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;

      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i , i) += p( i, i)*val ;

      } else if ( not ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i , j) += p( j, j)*val ;
        G( j , i) += p( j, j)*val ;
        G( j , j) += (p( j, i) + p( i, j))*val ;

      } else if ( not ieqj && i == k && j == l ){

        /* ( i j| i j) */

        G( i, j) += (p( j, i) + p( i, j))*val ;
        G( j, i) += (p( i, j) + p( j, i))*val ;

      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, i) += p( k, k)*val ;
        G( k, k) += p( i, i)*val ;

      } else if ( ieqj && not keql && i == k ){

        /* ( i i| i l) */

        G( i, i) += (p( i, l) + p( l, i))*val ;
        G( l, i) += p( i, i)*val ;
        G( i, l) += p( i, i)*val ;

      } else if ( not ieqj && not keql && j == l && i != k){

        /* ( i j| k j) */

        G( i, j) += (p( k, j) + p( j, k))*val ;
        G( j, i) += (p( k, j) + p( j, k))*val ;
        G( k, j) += (p( i, j) + p( j, i))*val ;
        G( j, k) += (p( i, j) + p( j, i))*val ;

      } else if ( not ieqj && keql && i != k && j != k) {

        /* ( i j| k k) */

        G( i, j) += p( k, k)*val ;
        G( j, i) += p( k, k)*val ;
        G( k, k) += (p( j, i) + p( i, j))*val ;

      } else if ( ieqj && not keql && i != k && i != l) {

        /* ( i i| k l) */

        G( i, i) += (p( k, l) + p( l, k))*val ;
        G( k, l) += p( i, i)*val ;
        G( l, k) += p( i, i)*val ;

      } else if ( not ieqj && not keql && j == k && i != l ) {

        /* ( i j| j l) */

        G( i, j)+= (p( j, l) + p( l, j))*val ;
        G( j, i)+= (p( j, l) + p( l, j))*val ;
        G( j, l)+= (p( i, j) + p( j, i))*val ;
        G( l, j)+= (p( i, j) + p( j, i))*val ;

      } else if ( not ieqj && not keql && i == k && j != l) {

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

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXcf> p, 
    Eigen::Ref<Eigen::MatrixXcf> G, const int nbasis) {

    /* Contract exchange integrals for complex matrices
    Given a spin density block contract the exchange terms into G */

    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    float val ;
    bool ieqk ;
    bool jeql ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() - 1 ;
      j = intarr[t].r_j() - 1 ;
      k = intarr[t].r_k() - 1 ;
      l = intarr[t].r_l() - 1 ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
      
      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i , i) += -p( i, i)*val ;

      } else if ( not ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i , j) += -p( j, j)*val ;
        G( j , i) += -p( j, j)*val ;
        G( j , j) += -(p( j, i) + p( i, j))*val ;

      } else if ( not ieqj && i == k && j == l ){

        /* ( i j| i J) */

        G( i, j) += -p( i, j)*val ;
        G( j, i) += -p( j, i)*val ;
        G( i, i) += -p( j, j)*val ;
        G( j, j) += -p( i, i)*val ;

      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, k) += -p( k, i)*val ;
        G( k, i) += -p( i, k)*val ;

      } else if ( ieqj && not keql && i == k ){

        /* ( i i| i l) */

        G( i, i) += -(p( i, l) + p( l, i))*val ;
        G( l, i) += -p( i, i)*val ;
        G( i, l) += -p( i, i)*val ;

      } else if ( not ieqj && not keql && j == l && i != k){

        /* ( i j| k j) */

        G( i, j) += -p( k, j)*val ;
        G( j, i) += -p( j, k)*val ;
        G( k, j) += -p( i, j)*val ;
        G( j, k) += -p( j, i)*val ;
        G( i, k) += -p( j, j)*val ;
        G( k, i) += -p( j, j)*val ;
        G( j, j) += -(p( i, k) + p( k, i))*val ;

      } else if ( not ieqj && keql && i != k ) {

        /* ( i j| k k) */

        G( i, k) += -p( k, j)*val ;
        G( j, k) += -p( k, i)*val ;
        G( k, j) += -p( i, k)*val ;
        G( k, i) += -p( j, k)*val ;

      } else if ( ieqj && not keql && i != k && i != l) {

        /* ( i i| k l) */

        G( i, l) += -p( k, i)*val ;
        G( i, k) += -p( l, i)*val ;
        G( k, i) += -p( i, l)*val ;
        G( l, i) += -p( i, k)*val ;

      } else if ( not ieqj && not keql && j == k ) {

        /* ( i j| j l) */

        G( i, l) += -p( j, j)*val ;
        G( l, i) += -p( j, j)*val ;
        G( l, j) += -p( i, j)*val ;
        G( j, l) += -p( j, i)*val ;
        G( i, j) += -p( l, j)*val ;
        G( j, i) += -p( j, l)*val ;
        G( j, j) += -(p( i, l) + p( l, i))*val ;

      } else if ( not ieqj && not keql && i == k && j != l) {

        /* ( i j| i l) */

        G( i, l) += -p( i, j)*val ;
        G( l, i) += -p( j, i)*val ;
        G( l, j) += -p( i, i)*val ;
        G( j, l) += -p( i, i)*val ;
        G( i, j) += -p( i, l)*val ;
        G( j, i) += -p( l, i)*val ;
        G( i, i) += -(p( j, l) + p( l, j))*val ;

      } else {

        /* ( i j| k l) */

        G( i, l)+= -p( k, j)*val ;
        G( j, l)+= -p( k, i)*val ;
        G( k, j)+= -p( i, l)*val ;
        G( j, k)+= -p( l, i)*val ;
        G( i, k)+= -p( l, j)*val ;
        G( l, j)+= -p( i, k)*val ;
        G( k, i)+= -p( j, l)*val ;
        G( l, i)+= -p( j, k)*val ;

      }

      }

    return 0 ;

  } ;

  int coulblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> G, const int nbasis) {

    /* Given a spin density block contract the coulomb terms into G 
     Contract coulomb integrals for real matrices */

    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    float val ;
    bool ieqk ;
    bool jeql ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() - 1 ;
      j = intarr[t].r_j() - 1 ;
      k = intarr[t].r_k() - 1 ;
      l = intarr[t].r_l() - 1 ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
      
      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i , i) += p( i, i)*val ;

      } else if ( not ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i , j) += p( j, j)*val ;
        G( j , i) += p( j, j)*val ;
        G( j , j) += (p( j, i) + p( i, j))*val ;

      } else if ( not ieqj && i == k && j == l ){

        /* ( i j| i J) */

        G( i, j) += (p( j, i) + p( i, j))*val ;
        G( j, i) += (p( i, j) + p( j, i))*val ;

      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, i) += p( k, k)*val ;
        G( k, k) += p( i, i)*val ;

      } else if ( ieqj && not keql && i == k ){

        /* ( i i| i l) */

        G( i, i) += (p( i, l) + p( l, i))*val ;
        G( l, i) += p( i, i)*val ;
        G( i, l) += p( i, i)*val ;

      } else if ( not ieqj && not keql && j == l && i != k){

        /* ( i j| k j) */

        G( i, j) += (p( k, j) + p( j, k))*val ;
        G( j, i) += (p( k, j) + p( j, k))*val ;
        G( k, j) += (p( i, j) + p( j, i))*val ;
        G( j, k) += (p( i, j) + p( j, i))*val ;

      } else if ( not ieqj && keql && i != k && j != k) {

        /* ( i j| k k) */

        G( i, j) += p( k, k)*val ;
        G( j, i) += p( k, k)*val ;
        G( k, k) += (p( j, i) + p( i, j))*val ;

      } else if ( ieqj && not keql && i != k && i != l) {

        /* ( i i| k l) */

        G( i, i) += (p( k, l) + p( l, k))*val ;
        G( k, l) += p( i, i)*val ;
        G( l, k) += p( i, i)*val ;

      } else if ( not ieqj && not keql && j == k ) {

        /* ( i j| j l) */

        G( i, j)+= (p( j, l) + p( l, j))*val ;
        G( j, i)+= (p( j, l) + p( l, j))*val ;
        G( j, l)+= (p( i, j) + p( j, i))*val ;
        G( l, j)+= (p( i, j) + p( j, i))*val ;

      } else if ( not ieqj && not keql && i == k && j != l) {

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

  int exchblk( std::vector<tei>& intarr, const Eigen::Ref<Eigen::MatrixXf> p, 
    Eigen::Ref<Eigen::MatrixXf> G, const int nbasis) {

    /* Given a spin density block contract the exchange terms into G 
     Contract exchange integrals for real matrices */

    int ntt ;
    int n2ei ;
    int i ;
    int j ;
    int k ;
    int l ;
    float val ;
    bool ieqk ;
    bool jeql ;
    bool ieqj ;
    bool keql ;

    ntt = nbasis*(nbasis + 1)/2 ;
    n2ei = ntt*(ntt+1)/2 ;

    for(int t = 0; t < n2ei; t++ ) {
      /*  ( 1 1| 2 2)
 *        ( i j| k l) = val */
      i = intarr[t].r_i() - 1 ;
      j = intarr[t].r_j() - 1 ;
      k = intarr[t].r_k() - 1 ;
      l = intarr[t].r_l() - 1 ;
      val = intarr[t].r_v() ;
//   
      ieqj = ( ( i == j) ? true : false) ;
      keql = ( ( k == l) ? true : false) ;
      
      if ( ieqj && keql && i == k ){

        /* ( i i| i i) */

        G( i , i) += -p( i, i)*val ;

      } else if ( not ieqj && keql && j == k ){

        /* ( i j| j j) */

        G( i , j) += -p( j, j)*val ;
        G( j , i) += -p( j, j)*val ;
        G( j , j) += -(p( j, i) + p( i, j))*val ;

      } else if ( not ieqj && i == k && j == l ){

        /* ( i j| i J) */

        G( i, j) += -p( i, j)*val ;
        G( j, i) += -p( j, i)*val ;
        G( i, i) += -p( j, j)*val ;
        G( j, j) += -p( i, i)*val ;

      } else if ( ieqj && keql && i != k ){

        /* ( i i| k k) */

        G( i, k) += -p( k, i)*val ;
        G( k, i) += -p( i, k)*val ;

      } else if ( ieqj && not keql && i == k ){

        /* ( i i| i l) */

        G( i, i) += -(p( i, l) + p( l, i))*val ;
        G( l, i) += -p( i, i)*val ;
        G( i, l) += -p( i, i)*val ;

      } else if ( not ieqj && not keql && j == l && i != k){

        /* ( i j| k j) */

        G( i, j) += -p( k, j)*val ;
        G( j, i) += -p( j, k)*val ;
        G( k, j) += -p( i, j)*val ;
        G( j, k) += -p( j, i)*val ;
        G( i, k) += -p( j, j)*val ;
        G( k, i) += -p( j, j)*val ;
        G( j, j) += -(p( i, k) + p( k, i))*val ;

      } else if ( not ieqj && keql && i != k ) {

        /* ( i j| k k) */

        G( i, k) += -p( k, j)*val ;
        G( j, k) += -p( k, i)*val ;
        G( k, j) += -p( i, k)*val ;
        G( k, i) += -p( j, k)*val ;

      } else if ( ieqj && not keql && i != k && i != l) {

        /* ( i i| k l) */

        G( i, l) += -p( k, i)*val ;
        G( i, k) += -p( l, i)*val ;
        G( k, i) += -p( i, l)*val ;
        G( l, i) += -p( i, k)*val ;

      } else if ( not ieqj && not keql && j == k ) {

        /* ( i j| j l) */

        G( i, l) += -p( j, j)*val ;
        G( l, i) += -p( j, j)*val ;
        G( l, j) += -p( i, j)*val ;
        G( j, l) += -p( j, i)*val ;
        G( i, j) += -p( l, j)*val ;
        G( j, i) += -p( j, l)*val ;
        G( j, j) += -(p( i, l) + p( l, i))*val ;

      } else if ( not ieqj && not keql && i == k && j != l) {

        /* ( i j| i l) */

        G( i, l) += -p( i, j)*val ;
        G( l, i) += -p( j, i)*val ;
        G( l, j) += -p( i, i)*val ;
        G( j, l) += -p( i, i)*val ;
        G( i, j) += -p( i, l)*val ;
        G( j, i) += -p( l, i)*val ;
        G( i, i) += -(p( j, l) + p( l, j))*val ;

      } else {

        /* ( i j| k l) */

        G( i, l)+= -p( k, j)*val ;
        G( j, l)+= -p( k, i)*val ;
        G( k, j)+= -p( i, l)*val ;
        G( j, k)+= -p( l, i)*val ;
        G( i, k)+= -p( l, j)*val ;
        G( l, j)+= -p( i, k)*val ;
        G( k, i)+= -p( j, l)*val ;
        G( l, i)+= -p( j, k)*val ;

      }

      }

    return 0 ;

  } ;

/* Evaluate elements between Slater Determinants */

float tranden ( common& com, hfwfn& a, hfwfn& b, Eigen::Ref<Eigen::MatrixXf> dabmat){

/*
 * Given two Slater Determinants, return the transition density in dabmat and return
 * the overlap.
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
 * Things that can be improved:
 *   - Add logic to find the transition density between different types
 *   of wavefunctions.
 *   - Allow for evaluating matrix elements between determinants in the non-orthogonal
 *   basis.
 *
 * */

  int wfn_a ;
  int wfn_b ;
  int nocc ;
  int nalp ;
  int nbet ;
  int nbasis ;
  int lstrt ;
  float ovl ;
  Eigen::MatrixXf moa ;
  Eigen::MatrixXf mob ;
  Eigen::MatrixXf dkl ;
  Eigen::MatrixXf tmp ;
  wfn_a = a.get_wti() ;
  wfn_b = b.get_wti() ;
  nocc = com.nele() ;
  dkl.resize( nocc, nocc) ;
  dkl.setZero() ;

  if ( wfn_a == 1 && wfn_b == 1 ) {

    nalp = com.nalp() ;
    nbasis = com.nbas() ;
    moa.resize( nbasis, nbasis) ;
    mob.resize( nbasis, nbasis) ;
    a.get_mos( moa) ;
    b.get_mos( mob) ;

    /* Build D*/
    dkl.block( 0, 0, nalp, nalp) = moa.block( 0, 0, nbasis, nalp).adjoint()*mob.block( 0, 0, nbasis, nalp) ;
    dkl.block( nalp, nalp, nalp, nalp) = dkl.block( 0, 0, nalp, nalp) ;

  } else if ( wfn_a == 3 && wfn_b == 3 ) {

    nalp = com.nalp() ;
    nbet = com.nbet() ;
    nbasis = com.nbas() ;
    moa.resize( nbasis, 2*nbasis) ;
    mob.resize( nbasis, 2*nbasis) ;
    a.get_mos( moa) ;
    b.get_mos( mob) ;
    dkl.block( 0, 0, nalp, nalp) = moa.block( 0, 0, nbasis, nalp).adjoint()*mob.block( 0, 0, nbasis, nalp) ;
    dkl.block( nalp, nalp, nbet, nbet) = moa.block( 0, nbasis, nbasis, com.nbet()).adjoint()*mob.block( 0, nbasis, nbasis, com.nbet()) ;

  } else if ( wfn_a == 5 && wfn_b == 5 ) {

    nbasis = 2*com.nbas() ;
    moa.resize( nbasis, nbasis) ;
    mob.resize( nbasis, nbasis) ;
    a.get_mos( moa) ;
    b.get_mos( mob) ;
    dkl = moa.block( 0, 0, nbasis, nocc).adjoint()*mob.block( 0, 0, nbasis, nocc) ;

  }

  /* Since all three will return the same size transition density handle the rest out here
   *
   * Build D(k | l)
   *
   * */

  ovl = dkl.determinant() ;

  /* Build the adjugate */
  tmp.resize( nocc, nocc) ;
  tmp = ovl*dkl.inverse() ;

  dkl = tmp ;
  tmp.resize( 0, 0) ;

  if ( wfn_a == 1 && wfn_b == 1 ) {
    dabmat.block( 0, 0, nbasis, nbasis) = mob.block( 0, 0, nbasis, nalp)*dkl.block( 0, 0, nalp, nalp)*moa.block( 0, 0, nbasis, nalp).adjoint() ;
    dabmat.block( nbasis, nbasis, nbasis, nbasis) = dabmat.block( 0, 0, nbasis, nbasis) ;
  } else if ( wfn_a == 3 && wfn_b == 3 ) {
    dabmat.block( 0, 0, nbasis, nbasis) = mob.block( 0, 0, nbasis, nalp)*dkl.block( 0, 0, nalp, nalp)*moa.block( 0, 0, nbasis, nalp).adjoint() ;
    dabmat.block( nbasis, nbasis, nbasis, nbasis) = mob.block( 0, nbasis, nbasis, nbet)*dkl.block( nalp, nalp, nbet, nbet)*moa.block( 0, nbasis, nbasis, nbet).adjoint() ;
  } else if ( wfn_a == 5 && wfn_b == 5 ) {
    dabmat = mob.block( 0, 0, nbasis, nocc)*dkl*moa.block( 0, 0, nbasis, nocc).adjoint() ;
  }

  mob.resize( 0, 0) ;
  moa.resize( 0, 0) ;
  dkl.resize( 0, 0) ;

  return ovl ;

} ;

std::complex<float> tranden ( common& com, hfwfn& a, hfwfn& b, Eigen::Ref<Eigen::MatrixXcf> dabmat) {

/* 
 *  This is an overloaded complex version of nointm
 * */

  int wfn_a ;
  int wfn_b ;
  int nocc ;
  int nalp ;
  int nbet ;
  int nbasis ;
  int lstrt ;
  cf ovl ;
  Eigen::MatrixXcf moa ;
  Eigen::MatrixXcf mob ;
  Eigen::MatrixXcf dkl ;
  Eigen::MatrixXcf tmp ;
  wfn_a = a.get_wti() ;
  wfn_b = b.get_wti() ;
  nocc = com.nele() ;
  dkl.resize( nocc, nocc) ;
  dkl.setZero() ;

  if ( wfn_a == 2 && wfn_b == 2 ) {

    nalp = com.nalp() ;
    nbasis = com.nbas() ;
    moa.resize( nbasis, nbasis) ;
    mob.resize( nbasis, nbasis) ;
    a.get_mos( moa) ;
    b.get_mos( mob) ;

    /* Build D*/
    dkl.block( 0, 0, nalp, nalp) = moa.block( 0, 0, nbasis, nalp).adjoint()*mob.block( 0, 0, nbasis, nalp) ;
    dkl.block( nalp, nalp, nalp, nalp) = dkl.block( 0, 0, nalp, nalp) ;

  } else if ( wfn_a == 4 && wfn_b == 4 ) {

    nalp = com.nalp() ;
    nbet = com.nbet() ;
    nbasis = com.nbas() ;
    moa.resize( nbasis, 2*nbasis) ;
    mob.resize( nbasis, 2*nbasis) ;
    a.get_mos( moa) ;
    b.get_mos( mob) ;
    dkl.block( 0, 0, nalp, nalp) = moa.block( 0, 0, nbasis, nalp).adjoint()*mob.block( 0, 0, nbasis, nalp) ;
    dkl.block( nalp, nalp, nbet, nbet) = moa.block( 0, nbasis, nbasis, com.nbet()).adjoint()*mob.block( 0, nbasis, nbasis, com.nbet()) ;

  } else if ( wfn_a == 6 && wfn_b == 6 ) {

    nbasis = 2*com.nbas() ;
    moa.resize( nbasis, nbasis) ;
    mob.resize( nbasis, nbasis) ;
    a.get_mos( moa) ;
    b.get_mos( mob) ;
    dkl = moa.block( 0, 0, nbasis, nocc).adjoint()*mob.block( 0, 0, nbasis, nocc) ;

  }

  /* Since all three will return the same size transition density handle the rest out here
   *
   * Build D(k | l) 
   *
   * */

  ovl = dkl.determinant() ;

  /* Build the adjugate */
  tmp.resize( nocc, nocc) ;
  tmp = ovl*dkl.inverse() ;

  dkl = tmp ;
  tmp.resize( 0, 0) ;

  if ( wfn_a == 2 && wfn_b == 2 ) {
    dabmat.block( 0, 0, nbasis, nbasis) = mob.block( 0, 0, nbasis, nalp)*dkl.block( 0, 0, nalp, nalp)*moa.block( 0, 0, nbasis, nalp).adjoint() ;
    dabmat.block( nbasis, nbasis, nbasis, nbasis) = dabmat.block( 0, 0, nbasis, nbasis) ;
  } else if ( wfn_a == 4 && wfn_b == 4 ) {
    dabmat.block( 0, 0, nbasis, nbasis) = mob.block( 0, 0, nbasis, nalp)*dkl.block( 0, 0, nalp, nalp)*moa.block( 0, 0, nbasis, nalp).adjoint() ;
    dabmat.block( nbasis, nbasis, nbasis, nbasis) = mob.block( 0, nbasis, nbasis, nbet)*dkl.block( nalp, nalp, nbet, nbet)*moa.block( 0, nbasis, nbasis, nbet).adjoint() ;
  } else if ( wfn_a == 6 && wfn_b == 6 ) {
    dabmat = mob.block( 0, 0, nbasis, nocc)*dkl*moa.block( 0, 0, nbasis, nocc).adjoint() ;
  }

  mob.resize( 0, 0) ;
  moa.resize( 0, 0) ;
  dkl.resize( 0, 0) ;

  return ovl ;

} ;

float obop ( common& com, Eigen::Ref<Eigen::MatrixXf> ouv, hfwfn& a, hfwfn& b) {

/*
 * One body operator.  Return <A|O|B> = sum_{ l k}<k|O|l> = sum{ u v} <u|O|v>C_{vl}D(k|l)C_{ku}^{*}
 * Passed two Slater determinants and matrix elements in an orthogonal basis, return the evaluated operator. */

  int nocc ;
  float aob ;
  Eigen::MatrixXf pvu ;
  Eigen::MatrixXf omega ;

  /* Build pvu */
 
  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.setZero() ;

  tranden ( com, a, b, pvu) ;

  omega = ouv*pvu ;
  aob = omega.trace() ;

  pvu.resize( 0, 0) ;
  omega.resize( 0, 0) ;

  return aob ;

} ;

std::complex<float> obop ( common& com, Eigen::Ref<Eigen::MatrixXcf> ouv, hfwfn& a, hfwfn& b) {

/*
 * Overloaded for complex functions.
 */

  int nbas ;

  std::complex<float> aob ;
  Eigen::MatrixXcf pvu ;
  Eigen::MatrixXcf omega ;

  /* Build pvu */

  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.setZero() ;

  tranden ( com, a, b, pvu) ;

  omega = ouv*pvu ;
  aob = omega.trace() ;

  return aob ;

} ;

float fockop ( common& com, Eigen::Ref<Eigen::MatrixXf> h, std::vector<tei>& intarr, hfwfn& a,
               hfwfn& b, float& ovl) { 

/*
 * While the fock operator is one body, it requires that we contract the two electron
 * integrals with the density to build it.  This wraps up that procedure.  This routine assumes
 *
 *  - Everything has been put into an orthogonal ao basis.  
 *  - Regardless of the flavor of hf wfn, this treats everything as ghf meaning 2*nbas dimensions.
 *  
 * */

  int nocc ;
  float aob ;
  Eigen::MatrixXf pvu ;
  Eigen::MatrixXf omega ;
  Eigen::MatrixXf f ;
  Eigen::MatrixXf g ;

  /* Build pvu */

  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
  f.resize( 2*com.nbas(), 2*com.nbas()) ;
  g.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.setZero() ;
  f.setZero() ;

  ovl = tranden ( com, a, b, pvu) ;

  ctr2eg( intarr, pvu, g, com.nbas()) ;

  f =  h + g ;
  g = h + f ;
  omega = g*pvu ;
  aob = 0.5*omega.trace() ;

  g.resize( 0, 0) ;
  f.resize( 0, 0) ;
  pvu.resize( 0, 0) ;
  omega.resize( 0, 0) ;

  return aob ;

} ;

std::complex<float> fockop ( common& com, Eigen::Ref<Eigen::MatrixXcf> h, std::vector<tei>& intarr, hfwfn& a, 
                             hfwfn& b, cf& ovl) {

/*
 * While the fock operator is one body, it requires that we contract the two electron
 * integrals with the density to build it.  This wraps up that procedure.  This routine assumes
 *
 *  - Everything has been put into an orthogonal ao basis.
 *  - Regardless of the flavor of hf wfn, this treats everything as ghf meaning 2*nbas dimensions.
 *
 * */

  int nocc ;
  cf aob ;
  cf pt5 = cf (0.5,0.0) ;
  Eigen::MatrixXcf pvu ;
  Eigen::MatrixXcf omega ;
  Eigen::MatrixXcf f ;
  Eigen::MatrixXcf fop ;
  Eigen::MatrixXcf g ;
  Eigen::MatrixXcf mos ;

  /* Build pvu */

  omega.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.resize( 2*com.nbas(), 2*com.nbas()) ;
  f.resize( 2*com.nbas(), 2*com.nbas()) ;
  g.resize( 2*com.nbas(), 2*com.nbas()) ;
  mos.resize( 2*com.nbas(), 2*com.nbas()) ;
  fop.resize( 2*com.nbas(), 2*com.nbas()) ;
  pvu.setZero() ;
  f.setZero() ;
  fop.setZero() ;
  a.get_mos( mos) ;

  ovl = tranden ( com, a, b, pvu) ;

  ctr2eg( intarr, pvu, g, com.nbas()) ;

  f =  h + g ;
  fop = mos.adjoint()*f*mos ;
  std::cout << " Fock matrix " << std::endl << fop << std::endl ;
  g = h + f ;
  omega = g*pvu ;
  aob = pt5*omega.trace() ;

  g.resize( 0, 0) ;
  f.resize( 0, 0) ;
  pvu.resize( 0, 0) ;
  omega.resize( 0, 0) ;

  return aob ;

} ;

