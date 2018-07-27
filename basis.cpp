#include "basis.h"
#include <string>
#include <vector>

 std::vector<libint2::Shell> build_basis ( std::string& bas_name, Eigen::Ref<Eigen::VectorXi> AtN, Eigen::Ref<Eigen::MatrixXd> coord) {
   /* Load a basis set */
   std::string b_sto3g = "sto-3g" ;

   if ( bas_name == b_sto3g ){
     sh = load_sto3g( AtN, coord) ;
   }

   return sh ;

 } ;

 std::vector<libint2::Shell> load_sto3g ( Eigen::VectorXi AtN, Eigen::MatrixXd c) {

  std::vector<libint2::Shell> shells;

  for( int a=0; a<AtN.size(); ++a) {

    // STO-3G basis set
    // cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
    //       doi: 10.1063/1.1672392
    // obtained from https://bse.pnl.gov/bse/portal
    switch (AtN(a)) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
              { // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}   // origin coordinates
            }
        );
        break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
              {71.616837000, 13.045096000, 3.530512200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        break;

      case 7: // Z=7: nitrogen
        shells.push_back(
            {
              {99.106169000, 18.052312000, 4.885660200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        shells.push_back(
            {
              {3.780455900, 0.878496600, 0.285714400},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        shells.push_back(
            {
          {3.780455900, 0.878496600, 0.285714400},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {130.709320000, 23.808861000, 6.443608300},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{c( a, 0), c( a, 1), c( a, 2)}}
            }
        );
        break;

      default:
        throw "do not know STO-3G basis for this Z";
    }

  }

  return shells;

 } ;
