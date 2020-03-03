/*
#include "common.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include "tei.h"
#include <vector>
*/
#include <string>

#ifndef QTZIO_H
#define QTZIO_H

namespace qtzio {

void parse_arguments(int c, char *v[]);
int read_input(const std::string& i);

}
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
