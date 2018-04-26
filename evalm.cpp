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

/*   */
