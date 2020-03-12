#include <Eigen/Dense>

namespace hamiltonian {

  class molecular {
    private:
      Eigen::VectorXd atom;
      Eigen::MatrixXd coord;
  
    public:
      molecular(void){};
      molecular(Eigen::VectorXd a, Eigen::MatrixXd c){
        atom = a;
        coord = c;
        };
      double nnrep(void);
  
  };
}
