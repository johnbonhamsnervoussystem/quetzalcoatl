#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <string>

#ifndef WFN_H
#define WFN_H

template<typename s, int r, int c>
class wfn {

  public :
  /* Track the type of wavefunction */
    int wfntyp ;
    double e_scf ;
  /* Store the wavefunction */
    Eigen::Matrix< s, r, c> moc ;
    Eigen::VectorXd eig ;
  /* Name for the binary file */
    std::string wfnfile ;
    wfn ( ): wfnfile("quetz.bqtz"){ }
    wfn ( std::string iofile): wfnfile(iofile){ }

} ;

template<typename s, int r, int c>
void save_wfn( wfn< s, r, c>& w, int cntl = 0) ;

template<typename s, int r, int c>
void load_wfn( wfn< s, r, c>& w, int cntl = 0) ;

#endif
