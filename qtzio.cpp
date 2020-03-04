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
#include "periodic_table.h"
#include "qtzio.h"
#include <string>
#include <vector>


/* Routines to swtich between string options */

namespace qtzio {

void parse_arguments(int argc, char *argv[]){
  /*
   * Ensure that the command line arguments are valid and exit if they are incorrect.
   */
  if (argc != 2){
    std::cout << "Expected 1 argument but received " << (argc - 1) << std::endl;
    std::exit(EXIT_FAILURE);
    }

  std::string extension = "";
  std::string filename(*(argv + 1));
  auto itr = filename.begin();

  for (auto i = (filename.length() - 5); i < filename.length(); i++){
    extension += *(itr + i);
    }

  if (extension != ".json"){
    std::cout << "Expected a .json file but received " << extension << std::endl;
    std::exit(EXIT_FAILURE);
    }

  std::cout << "Input file: " << filename << std::endl;

  qtzio::read_input(filename);

  return;

  }

int read_input(const std::string& inpfile){
  /*
   * Read the input file and save the information into common.
   */
  int natoms;
  std::vector<double> atoms;
  std::vector<std::vector<double>> coordinates;
  Json::Value root;
  std::ifstream ifs;
  std::string data;
  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;

  ifs.open(inpfile, std::ifstream::in);
  if (!parseFromStream(builder, ifs, &root, &errs)) {
    std::cout << errs << std::endl;
    return EXIT_FAILURE;
    }
  natoms = root["geometry"]["atoms"].size();
  atoms.reserve(natoms);
  coordinates.reserve(natoms);

  for (Json::Value::iterator itr = root["geometry"]["atoms"].begin(); itr != root["geometry"]["atoms"].end(); itr++){
    atoms.push_back(static_cast<double>(periodic_table::symbol_conversion[boost::algorithm::to_lower_copy((*itr).asString())]));
    }
//  std::cout.precision(2);
//  for (std::vector<double>::iterator itr = atoms.begin(); itr != atoms.end(); itr++){
//    std::cout << std::fixed << *itr << std::endl;
//    }

  for (Json::Value::iterator itr = root["geometry"]["coordinates"].begin(); itr != root["geometry"]["coordinates"].end(); itr++){
    std::vector<double> tmp;
    for (unsigned int itc = 0; itc < (*itr).size(); itc ++){
      tmp.push_back((*itr)[itc].asDouble());
      }
    coordinates.push_back(tmp);
    }

  qtzio::print_system(atoms, coordinates);

  return EXIT_SUCCESS;

};

void print_system(std::vector<double> a, std::vector<std::vector<double>> c){
  std::vector<double>::iterator itr;
  std::cout << " atoms ";
  std::cout << "    x       ";
  std::cout << "    y       ";
  std::cout << "    z       ";
  std::cout << std::endl;
  for (auto atom = 0; atom < a.size(); atom++){
    std::cout << " " << std::setw(2) << periodic_table::atomic_number_conversion[a[atom]];
    for (auto itr = c[atom].begin(); itr != c[atom].end(); itr++){
      std::cout << std::setw(11) << *itr << " ";
      }
    std::cout << std::endl;
    }
  return;
}

}



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
