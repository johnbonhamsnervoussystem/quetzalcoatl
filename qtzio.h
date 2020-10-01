/*
#include "common.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include "tei.h"
*/
#include <json.h>
#include <string>
#include <vector>

#include <libint2.hpp>
#include "qtzctl.h"

#ifndef QTZIO_H
#define QTZIO_H


class QtzInput {
  private:
    std::string inputfile;
    Json::Value root_input;
    void check_members(void);
    void parse_method();
    void parse_molecular_input(Json::Value m);
    void print_system(std::vector<double> a, std::vector<std::vector<double>> c);
  public:
    QtzInput(int argc, char *argv[]);
    void parse_input(void);
    QtzControl control(void);
    std::vector<libint2::Atom> atoms(void);
  };
//int read_input(const std::string& i);
//void parse_molecular_input(const Json::Value r);
//void print_system(std::vector<double> a, std::vector<std::vector<double>> c);

/*
bool open_text( std::ofstream& F_OUT, int cntl, const std::string& f = "qtztemp") ;

bool open_binary( std::ofstream& F_OUT, int cntl) ;

bool open_binary( std::ifstream& F_IN) ;

void strip_lower( std::string& s) ;

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ofstream& F_OUT) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::ifstream& F_IN) ;

template <class matrix>
void print_mat( const matrix& o, std::string h = "") ;
*/
#endif
