/* Routines for Quetz I/O and printing
#include "common.h"
#include <complex>
#include "constants.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include "qtzcntrl.h"
#include <sys/stat.h>
#include "tei.h"
#include "time_dbg.h"
#include <unordered_map>
*/
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <json.h>
#include <libint2.hpp>
#include "periodic_table.h"
#include "qtzctl.h"
#include "qtzio.h"
#include <string>
#include <vector>


QtzInput::QtzInput(int argc, char *argv[]) {
  if (argc != 2){
    std::cout << "Expected 1 argument but received " << (argc - 1) << std::endl;
    std::exit(EXIT_FAILURE);
    }
  
  inputfile = *(argv + 1);
  std::string::iterator itr = inputfile.begin();
  
  std::string extension = "";
  for (auto i = (inputfile.length() - 5); i < inputfile.length(); i++){
    extension += *(itr + i);
    }
  
  if (extension != ".json"){
    std::cout << "Expected a .json file but received " << extension << std::endl;
    std::exit(EXIT_FAILURE);
    }
  
  std::cout << "input file: " << inputfile << std::endl;
  };

QtzControl QtzInput::control(void){
  QtzControl qtz_control;
  qtz_control.directive = "wavefunction";
  return qtz_control;
  }

std::vector<libint2::Atom> QtzInput::atoms(void){
    int natoms = root_input["geometry"]["atoms"].size();
    Json::Value::iterator itr;
    std::vector<libint2::Atom> atoms;
    std::vector<libint2::Atom>::iterator ita;
    for (auto itn = 0; itn < natoms; itn++) {atoms.push_back(libint2::Atom());}

    ita = atoms.begin();
    for (itr = root_input["geometry"]["atoms"].begin(); itr != root_input["geometry"]["atoms"].end(); itr++){
      ita->atomic_number = static_cast<double>(periodic_table::symbol_conversion[boost::algorithm::to_lower_copy((*itr).asString())]);
      ita++;
      }
  
    ita = atoms.begin();
    for (itr = root_input["geometry"]["coordinates"].begin(); itr != root_input["geometry"]["coordinates"].end(); itr++){
       std::cout << *itr << std::endl;
/*
      ita->x = static_cast<double>(itr[0]);
      ita->y = static_cast<double>(itr[1]);
      ita->z = static_cast<double>(itr[2]);
      ita++;
*/
      }

    return atoms;

  };

void QtzInput::parse_input(void){
  Json::Value root_input;
  Json::CharReaderBuilder builder;
  Json::Value::iterator itr;
  JSONCPP_STRING errs;
  std::ifstream ifs;
  
  ifs.open(inputfile, std::ifstream::in);
  if (!parseFromStream(builder, ifs, &root_input, &errs)) {
    std::cout << errs << std::endl;
    std::exit(EXIT_FAILURE);
    }

    check_members(root_input);
    parse_method(root_input["method"]);

  if (root_input["hamiltonian"] == "molecular") {
    parse_molecular_input(root_input["geometry"]);
    };

  return;

  };

void QtzInput::check_members(Json::Value input_json){
  /*
   * Check that the json input has all the required members.
   * */
  bool error = false;
  std::vector<std::string>::iterator itr;
  std::vector<std::string> molecular_requirements{"geometry", "basis-set"};
  if (input_json["hamiltonian"] == "molecular") {
    for (itr = molecular_requirements.begin(); itr != molecular_requirements.end(); itr++){
      if (not input_json.isMember(*itr)){
        error = true;
        std::cout << *itr << " is missing" << std::endl;
        }
      }
    }
  if (error){
    std::exit(EXIT_FAILURE);
    }

  return;

  }

void QtzInput::parse_method(Json::Value method){
  if (method == "rhf") {
    std::cout << "Restricted Hartree-Fock" << std::endl;
  } else {
    std::cout << "Unrecognized method" << std::endl;
    std::exit(EXIT_FAILURE);
    }

  return;

  }

void QtzInput::parse_molecular_input(Json::Value molecule){
    int natoms = molecule["atoms"].size();
    Json::Value::iterator itr;
    std::vector<double> atoms;
    atoms.reserve(natoms);
    std::vector<std::vector<double>> coordinates;
    coordinates.reserve(natoms);
  
    for (itr = molecule["atoms"].begin(); itr != molecule["atoms"].end(); itr++){
      atoms.push_back(static_cast<double>(periodic_table::symbol_conversion[boost::algorithm::to_lower_copy((*itr).asString())]));
      }
  
    for (itr = molecule["coordinates"].begin(); itr != molecule["coordinates"].end(); itr++){
      std::vector<double> tmp;
      for (unsigned int itc = 0; itc < (*itr).size(); itc ++){
        tmp.push_back((*itr)[itc].asDouble());
        }
      coordinates.push_back(tmp);
      }

    print_system(atoms, coordinates);

    return;

  };

  void QtzInput::print_system(std::vector<double> a, std::vector<std::vector<double>> c){
    std::vector<double>::iterator itr;
    std::cout << " atoms ";
    std::cout << "   _x_      ";
    std::cout << "   _y_      ";
    std::cout << "   _z_      ";
    std::cout << std::endl;
    for (auto atom = 0; atom < a.size(); atom++){
      std::cout << " " << std::setw(2) << periodic_table::atomic_number_conversion[a[atom]];
      for (auto itr = c[atom].begin(); itr != c[atom].end(); itr++){
        std::cout << std::setw(11) << *itr << " ";
        }
      std::cout << std::endl;
      }
    return;
  };

//
//bool open_text( std::ofstream& F_OUT, int cntl, const std::string& filename) {
///*
//  cntl :
//    0 - (default) write to a new file or overwrite one that already exists
//    1 - Check if a file exists and append.  If it does not exist create it.
//*/
//  const std::string wfnIO = filename + ".txt" ;
//  bool ioerr = false ;
//  struct stat buf ;
//
//  if ( cntl == 0 ){
//    try {
//      F_OUT.open( wfnIO, std::ofstream::out) ;
//    } catch ( std::fstream::failure& e) {
//      std::cout << e.what() << std::endl ;
//      ioerr = true ;
//    } catch (...) {
//      std::cout << " Default excepetion inside write_to_bin " << std::endl ;
//      ioerr = true ;
//      }
//  } else if ( cntl == 1 ){
//    /* Check if the file exists */
//    if ( stat(wfnIO.c_str(), &buf) == 0) {
//      try {
//        F_OUT.open( wfnIO, std::ofstream::app | std::ofstream::out) ;
//      } catch ( std::fstream::failure& e) {
//        std::cout << e.what() << std::endl ;
//        ioerr = true ;
//      } catch (...) {
//        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
//        ioerr = true ;
//        }
//    } else {
//      try {
//        F_OUT.open( wfnIO, std::ofstream::out) ;
//      } catch ( std::fstream::failure& e) {
//        std::cout << e.what() << std::endl ;
//        ioerr = true ;
//      } catch (...) {
//        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
//        ioerr = true ;
//        }
//    }
//  } else {
//    qtzcntrl::shutdown( "Unrecognized Option (open_binary)") ;
//    }
//
//  return ioerr ;
//
//} ;
//
//bool open_binary( std::ofstream& F_OUT, int cntl) {
///*
//  cntl :
//    0 - (default) write to a new file or overwrite one that already exists
//    1 - Check if a file exists.  If not create it.
//*/
//  const std::string wfnIO = "qtz.wfn.bin" ;
//  bool ioerr = false ;
//  struct stat buf ;
//
//  if ( cntl == 0 ){
//    try {
//      F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
//    } catch ( std::fstream::failure& e) {
//      std::cout << e.what() << std::endl ;
//      ioerr = true ;
//    } catch (...) {
//      std::cout << " Default excepetion inside write_to_bin " << std::endl ;
//      ioerr = true ;
//      }
//  } else if ( cntl == 1 ){
//    /* Check if the file exists */
//    if ( stat(wfnIO.c_str(), &buf) == 0) {
//      try {
//        F_OUT.open( wfnIO, std::ofstream::app | std::ofstream::binary | std::ofstream::out) ;
//      } catch ( std::fstream::failure& e) {
//        std::cout << e.what() << std::endl ;
//        ioerr = true ;
//      } catch (...) {
//        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
//        ioerr = true ;
//        }
//    } else {
//      try {
//        F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
//      } catch ( std::fstream::failure& e) {
//        std::cout << e.what() << std::endl ;
//        ioerr = true ;
//      } catch (...) {
//        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
//        ioerr = true ;
//        }
//    }
//  } else {
//    qtzcntrl::shutdown( "Unrecognized Option (open_binary)") ;
//    }
//
//  return ioerr ;
//
//} ;
//
//bool open_binary( std::ifstream& F_IN) {
///*
//  (default) Check if a file exists and if it does open it
//*/
//  const std::string wfnIO = "qtz.wfn.bin" ;
//  bool ioerr = false ;
//  struct stat buf ;
//
//  /* Check if the file exists */
//  if ( stat(wfnIO.c_str(), &buf) == 0) {
//    try {
//      F_IN.open( wfnIO, std::ofstream::binary | std::ofstream::in) ;
//    } catch ( std::fstream::failure& e) {
//      std::cout << e.what() << std::endl ;
//      ioerr = true ;
//    } catch (...) {
//      std::cout << " Default exception inside write_to_bin " << std::endl ;
//      ioerr = true ;
//      }
//  } else {
//    qtzcntrl::shutdown( "FIle does not exist (open_binary)") ;
//    }
//
//  return ioerr ;
//
//} ;
//
//void strip_lower( std::string& s) {
///*
//  Wrap together whitespace trimming and lower casing for input file parsing
//*/
//  std::string::iterator end_pos ;
//  end_pos = std::remove( s.begin(), s.end(), ' ') ;
//  s.erase( end_pos, s.end()) ;
//  std::transform(s.begin(), s.end(), s.begin(), ::tolower) ;
//
//  return ;
//
//} ;
//
//template <class matrix> 
//void write_eigen_bin (const matrix& m, std::ofstream& F_OUT) {
//  typename matrix::Index rows=m.rows(), cols=m.cols();
//
//  F_OUT.write( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;
//
//  return ;
//
//}
//
//template void write_eigen_bin(const Eigen::VectorXd&, std::ofstream&) ;
//template void write_eigen_bin(const Eigen::VectorXcd&, std::ofstream&) ;
//template void write_eigen_bin(const Eigen::MatrixXd&, std::ofstream&) ;
//template void write_eigen_bin(const Eigen::MatrixXcd&, std::ofstream&) ;
//
//template <class matrix> 
//void read_eigen_bin (const matrix& m, std::ifstream& F_IN) {
//  typename matrix::Index rows=m.rows(), cols=m.cols();
//
//  F_IN.read( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;
//
//  return ;
//
//}
//
//template void read_eigen_bin(const Eigen::VectorXd&, std::ifstream&) ;
//template void read_eigen_bin(const Eigen::VectorXcd&, std::ifstream&) ;
//template void read_eigen_bin(const Eigen::MatrixXd&, std::ifstream&) ;
//template void read_eigen_bin(const Eigen::MatrixXcd&, std::ifstream&) ;
//
//template <class matrix>
//void print_mat( const matrix& o, std::string h){
//  int i = 0 ;
//  typename matrix::Index cols = o.cols(), rows = o.rows() ;
//  Eigen::IOFormat matprt(5, 0, "  ", "\n", "| ", "|", "") ;
//
//  if ( ! h.empty()){
//    std::cout << h << std::endl ;
//    std::cout << "---- ---- ----" << std::endl ;
//    }
//
//  if ( cols <= 5){
//      std::cout << o.format( matprt) << std::endl ;
//  } else {
//    do {
//      std::cout << o.block( 0, i, rows, 5).format( matprt) << std::endl << std::endl ;
//      i += 5 ;
//      cols -= 5 ;
//      } while ( cols > 5) ;
//      std::cout << o.block( 0, i, rows, cols).format( matprt) << std::endl ;
//    }
//
//  return ;
//
//}
//
//template void print_mat(const Eigen::VectorXd& , std::string) ;
//template void print_mat(const Eigen::VectorXcd& , std::string) ;
//template void print_mat(const Eigen::MatrixXd& , std::string) ;
//template void print_mat(const Eigen::MatrixXcd& , std::string) ;
//template void print_mat(const Eigen::Ref<Eigen::VectorXd>& , std::string) ;
//template void print_mat(const Eigen::Ref<Eigen::VectorXcd>& , std::string) ;
//template void print_mat(const Eigen::Ref<Eigen::MatrixXd>& , std::string) ;
//template void print_mat(const Eigen::Ref<Eigen::MatrixXcd>& , std::string) ;
//
