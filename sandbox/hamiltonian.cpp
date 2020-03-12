#include "hamiltonian.h"
#include <iostream>
#include <Eigen/Dense>

#define d0 0.e0
#define d2 2.e0

namespace hamiltonian {

    double molecular::nnrep(void) {
      /*
        Get the nuclear-nuclear repulsion
      */
      int i, j, natm = atom.size();
      double cx = d0, cy = d0, cz = d0;
      double r2 = d0, n_rep = d0;

      for (i = 0; i < natm; i++){
        for (j = i+1; j < natm; j++){
          cx = coord(j, 0) - coord(i, 0);
          cy = coord(j, 1) - coord(i, 1);
          cz = coord(j, 2) - coord(i, 2);
          r2 = std::pow(cx, d2) + std::pow(cy, d2) + std::pow(cz, d2);
          n_rep += atom(i)*atom(j)/std::sqrt(r2);
          }
        }
    
      std::cout << " Nuclear repulsion is " << n_rep << std::endl;
    
      return n_rep;
    
      };

};
