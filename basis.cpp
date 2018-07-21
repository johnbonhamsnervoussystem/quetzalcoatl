#include "basis.h"
#include <libint2.h>
#include <string>
#include <vector>

 std::vector<libint2::Shell> build_basis ( std::string& bas_name, Eigen::Ref<Eigen::VectorXi> AtN, Eigen::Ref<Eigen::MatrixXd> coord) {
   /* Load a basis set */
   std::string b_sto3g = "sto-3g" ;
   std::vector<libint2::Shell> sh ;

   if ( bas_name == b_sto3g ){
     sh = load_sto3g( Eigen::VectorXi AtN, Eigen::MatrixXd coord) ;
   }

   return sh ;

 } ;

