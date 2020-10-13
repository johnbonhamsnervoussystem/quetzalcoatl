#ifndef PAIRING_H
#define PAIRING_H

template <class matrix>
class pairing: public nbodyint<matrix> {
  private :
    using nbodyint<matrix>::itype ;
    using nbodyint<matrix>::dim ;
    using nbodyint<matrix>::G ;
    double U ;
  public :
    pairing( double U, int i, int n) : nbodyint<matrix>::nbodyint( i, n){
      this->U = U ;
      } ;
    ~pairing(){} ;
    void contract( matrix& m) ;
    void contract( matrix& m, matrix& n) ;
} ;

#endif
