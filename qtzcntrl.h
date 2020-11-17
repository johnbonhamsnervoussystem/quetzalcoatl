#include <libint2.hpp>
#include <vector>
#include <string>

#ifndef QTZCNTRL_H
#define QTZCNTRL_H

class QtzControl {
  public:
    std::string basis;
    std::string directive;
    std::string hamiltonian;
    std::string method;
    std::vector<libint2::Atom> atoms;
    void shutdown(std::string s);
};

#endif
