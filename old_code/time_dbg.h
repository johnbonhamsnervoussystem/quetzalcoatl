#include <chrono>
#include <string>

#ifndef TIMER_H
#define TIMER_H

class time_dbg {
  private :
    std::string r_name ;
    std::chrono::system_clock::time_point start ;

  public :
    time_dbg ( std::string s) ;
    void end( void) ;

} ;
#endif
