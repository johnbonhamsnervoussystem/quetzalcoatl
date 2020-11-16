#include <Eigen/Dense>

namespace hamiltonian {

  class molecular {
    private:
      Eigen::VectorXd atom;
      Eigen::MatrixXd coord;
      Eigen::MatrixXd one_body;
      Eigen::MatrixXd potenial;
      Eigen::MatrixXd kinetic;
  
    public:
      molecular(void){};
      molecular(Eigen::VectorXd a, Eigen::MatrixXd c){
        atom = a;
        coordinates = c;
        };

      void coordinates(Eigen::MatrixXd c);
      void atomic_number(Eigen::VectorXd a);
      double nnrep(void);
  
  };
}
