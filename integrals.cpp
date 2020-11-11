#include <iostream>
#include <libint2.hpp>
#include "integrals.h"
#include <Eigen/Core>


void compute_overlap(libint2::Engine engine, libint2::BasisSet basis) {
  auto shell2bf = basis.shell2bf();
  const auto& buf_vec = engine.results();
  
  for(auto s1=0; s1!=basis.size(); ++s1) {
    for(auto s2=0; s2!=basis.size(); ++s2) {
  
      std::cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
      engine.compute(basis[s1], basis[s2]);
      std::cout << "done" << std::endl;
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;
  
      auto bf1 = shell2bf[s1];
      auto n1 = basis[s1].size();
      auto bf2 = shell2bf[s2];
      auto n2 = basis[s2].size();
  
      for(auto f1=0; f1!=n1; ++f1)
        for(auto f2=0; f2!=n2; ++f2)
          std::cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << std::endl;
    }
  }
  return;
};
