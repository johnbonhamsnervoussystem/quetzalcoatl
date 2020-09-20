#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>

#ifndef WFN_H
#define WFN_H

template<typename s, int r, int c>
class wfn {

  public :
  /* Store the wavefunction */
    Eigen::Matrix< s, r, c> moc ;

} ;

#endif
