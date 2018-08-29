#include <chrono>
#include <iostream>
#include <string>
#include "time_dbg.h"

  time_dbg::time_dbg ( std::string s){
    start = std::chrono::system_clock::now() ;
    r_name = s ;
    std::cout << "Entering " << s << std::endl ;
    } ;
 
  void time_dbg::end( void){
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now() ;
    std::chrono::duration<double> span = end - start ;
    std::cout << r_name << " --time-- " << span.count() << std::endl ;
    return ;
    } ;
