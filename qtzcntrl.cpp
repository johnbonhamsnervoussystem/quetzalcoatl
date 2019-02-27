#include <iostream>
#include "qtzcntrl.h"
#include <string>

void qtzcntrl::shutdown( std::string s) {
/*
  Basic shutdown function.
*/
    std::cout << s << std::endl ;
    std::cout << "Exiting Quetzalcoatl" << std::endl ;
    exit(EXIT_FAILURE) ;

    return ;

  } 

void qtzcntrl::end( std::string s) {
/*
  Super basic shutdown function.
*/
    std::cout << s << std::endl ;
    exit(EXIT_SUCCESS) ;

    return ;

  } 
