#include <iostream>
#include <libint2.hpp>
#include "integrals.h"
#include "qtzctl.h"
#include "qtzio.h"
#include <string>
#include <vector>
#include <Eigen/Core>

int main(int argc, char *argv[]) {
  
  QtzInput parser(argc, argv);
  parser.parse_input();
  QtzControl qtz_control = parser.control();
  std::vector<libint2::Atom> atoms = parser.parse_atoms();

  /*
   * Figure out what we need to do.  Let's start with solving for wavefunction
   * */
  if (qtz_control.directive == "wavefunction") {
    std::cout << "Solving for a wavefunction" << std::endl;
    libint2::initialize();
    libint2::BasisSet basis_set(parser.basis_set(), atoms);
    int nbasis = basis_set.nbf();
    std::cout << "Number of basis functions " << nbasis << std::endl;
    libint2::Engine s_engine(libint2::Operator::overlap, basis_set.max_nprim(), basis_set.max_l());


    compute_overlap(s_engine, basis_set);

    libint2::finalize();
    }

  return 0 ;

}

