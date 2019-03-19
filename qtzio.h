#include "common.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include "tei.h"
#include <vector>

#ifndef QTZIO_H
#define QTZIO_H

void read_input( common& c, const std::string& in_s) ;

bool open_text( std::ofstream& F_OUT, int cntl, const std::string& f = "qtztemp") ;

bool open_binary( std::ofstream& F_OUT, std::string& wfnIO, int cntl) ;

bool open_binary( std::ifstream& F_IN, std::string& wfnIO) ;

template <class matrix>
void print_mat( const matrix& o, std::string h = "") ;

void process_input_line( std::string& r, std::vector<std::string>& O) ;

void strip_lower( std::string& s) ;

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ofstream& F_OUT) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::ifstream& F_IN) ;

#endif
