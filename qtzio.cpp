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
  
  this->inputfile = *(argv + 1);
  std::string::iterator itr = this->inputfile.begin();
  
  std::string extension = "";
  for (auto i = (this->inputfile.length() - 5); i < this->inputfile.length(); i++){
    extension += *(itr + i);
    }
  
  if (extension != ".json"){
    std::cout << "Expected a .json file but received " << extension << std::endl;
    std::exit(EXIT_FAILURE);
    }
  
  std::cout << "input file: " << this->inputfile << std::endl;
  };


std::string QtzInput::basis_set(void) {
  Json::StreamWriterBuilder builder;
  std::string source = Json::writeString(builder, this->root_input["basis-set"]);
  source.erase(remove(source.begin(), source.end(), '\"'), source.end());
  return source;
};


std::string QtzInput::hamiltonian(void) {
  Json::StreamWriterBuilder builder;
  std::string source = Json::writeString(builder, this->root_input["basis-set"]);
  source.erase(remove(source.begin(), source.end(), '\"'), source.end());
  return source;
};


QtzControl QtzInput::control(void){
  QtzControl qtz_control;
  qtz_control.basis = this->basis_set();
  qtz_control.directive = "wavefunction";
  qtz_control.hamiltonian = this->hamiltonian();
  qtz_control.method = this->parse_method();
  qtz_control.atoms = this->parse_atoms();
  return qtz_control;
  }

std::vector<libint2::Atom> QtzInput::parse_atoms(void){
    Json::Value::iterator itr;
    std::vector<libint2::Atom> atoms;
    std::vector<libint2::Atom>::iterator ita;
    Json::Value atom_block = this->root_input["geometry"]["atoms"];
    Json::Value coordinate_block = this->root_input["geometry"]["coordinates"];
    int natoms = atom_block.size();
    for (auto itn = 0; itn < natoms; itn++) {atoms.push_back(libint2::Atom());}

    ita = atoms.begin();
    for (itr = atom_block.begin(); itr != atom_block.end(); itr++) {
      ita->atomic_number = static_cast<double>(periodic_table::symbol_conversion[boost::algorithm::to_lower_copy((*itr).asString())]);
      ita++;
      }
  
    ita = atoms.begin();
    for (itr = coordinate_block.begin(); itr != coordinate_block.end(); itr++) {
      for (unsigned int itc = 0; itc < (*itr).size(); itc++) {
	if (itc == 0) {ita->x = (*itr)[itc].asDouble();}
	if (itc == 1) {ita->y = (*itr)[itc].asDouble();}
	if (itc == 2) {ita->z = (*itr)[itc].asDouble();}
        }
      ita++;
      }

    return atoms;

  };

void QtzInput::parse_input(void){
  Json::CharReaderBuilder builder;
  Json::Value::iterator itr;
  JSONCPP_STRING errs;
  std::ifstream ifs;
  
  ifs.open(this->inputfile, std::ifstream::in);
  if (!parseFromStream(builder, ifs, &this->root_input, &errs)) {
    std::cout << errs << std::endl;
    std::exit(EXIT_FAILURE);
    }

  this->check_members();
  this->parse_method();

  if (root_input["hamiltonian"] == "molecular") {
    parse_molecular_input(this->root_input["geometry"]);
    };
  return;
  };

void QtzInput::check_members(void){
  /*
   * Check that the json input has all the required members.
   * */
  bool error = false;
  std::vector<std::string>::iterator itr;
  std::vector<std::string> molecular_requirements{"geometry", "basis-set"};
  if (this->root_input["hamiltonian"] == "molecular") {
    for (itr = molecular_requirements.begin(); itr != molecular_requirements.end(); itr++){
      if (not this->root_input.isMember(*itr)){
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

void QtzInput::parse_method(void){
  if (this->root_input["method"] == "rhf") {
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

int main(int argc, char* argv []) {

  QtzInput parser(argc, argv);
  parser.parse_input();
  QtzControl qtz_control = parser.control();
  return 0;

}
