#include <iostream>
#include <string>
#include "qtzscratch.h"
#include <Eigen/Core>


bool QtzScratch::is_member(std::string key) {
  return this->_doublemat.count(key) == 1;
}


void QtzScratch::new_dblmat(std::string key, int nbasis) {
  Eigen::MatrixXd* dmp = new Eigen::MatrixXd;
  dmp->resize(nbasis, nbasis);
  this->_doublemat.insert({key, dmp});
}


void QtzScratch::print(std::string key) {
  if (this->is_member(key)){
    std::cout << *(this->_doublemat[key]) << std::endl;
  } else {
    std::cout << "not a member" << std::endl;
    }
}


Eigen::MatrixXd* QtzScratch::dblmat(std::string key) {
  std::cout << "returning pointer to matrix" << std::endl;
  return this->_doublemat[key];
}


int main(int argc, char* argv[]){
  int nbasis = 2;
  QtzScratch a;
  a.new_dblmat("overlap", nbasis);
  if (a.is_member("overlap")) {
    std::cout << "overlap is member" << std::endl;
    }
  std::cout << "accessing pointer" << std::endl;
  Eigen::MatrixXd* s = a.dblmat("overlap");
  (*s)(0, 0) = 1.0;
  (*s)(1, 0) = 2.0;
  (*s)(0, 1) = 3.0;
  (*s)(1, 1) = 4.0;
  a.print("overlap");

  return 0;

  }
