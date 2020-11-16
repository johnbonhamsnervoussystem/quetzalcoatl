#include "tei.h"

  void tei::set(int a, int b, int c, int d, float integral ){
    i = a ;
    j = b ;
    k = c ;
    l = d ;
    val = integral ;
  }

  float tei::r_v ( void) { return val ;}
  int tei::r_i ( void) { return i ;}
  int tei::r_j ( void) { return j ;}
  int tei::r_k ( void) { return k ;}
  int tei::r_l ( void) { return l ;}

