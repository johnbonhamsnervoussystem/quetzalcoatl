#include <iostream>
#include "qtzcntrl.h"
#include <string>

void QtzControl::shutdown(std::string s) {
/*
  Super basic shutdown function.
*/
    std::cout << s << std::endl ;
    std::cout << "Exiting Quetzalcoatl" << std::endl ;
    exit(EXIT_FAILURE) ;

    return ;

  }
