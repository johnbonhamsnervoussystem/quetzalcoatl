/* Routines for Quetz I/O and printing */
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "common.h"
#include <complex>
#include "constants.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "qtzcntrl.h"
#include "qtzio.h"
#include <sys/stat.h>
#include <string>
#include "tei.h"
#include "time_dbg.h"
#include <unordered_map>
#include <vector>

/*

  Reading and parsing the input for traversal through the program may end
  up being the most complicated part of the whole program

  This needs to be clear and easy to modify which is currently is not. 
  Unfortunately this is likely the last thing to be cleaned.

*/

void read_input( common& com, const std::string& inpfile){
/* 
  Read the input file and save the information into common.
*/
  int i_junk, t_junk ;
  int at_count = 0 ;
  std::vector<double> tmp(3) ;
  double d_junk, d_calc ;
  size_t pos = 0 ;
  std::string line ;
  std::string s_junk, s_j1, s_j2 ;
  std::string delim = "," ;
  std::string::size_type first, last ;
  std::ifstream jobfile ;
  std::vector<double> atnum ;
  std::vector<std::vector<double>> t_c ;
  std::vector<std::string> OUT, OUTx ;
  std::vector<std::string>::iterator it ;
  std::unordered_map< std::string, int> prjcase ( { {"n", 1}, {"s", 2}, {"t", 3}, {"k", 4}} ) ;
  time_dbg read_input_time = time_dbg("read_input") ;

  jobfile.open( inpfile, std::ifstream::in ) ;

  while( getline( jobfile, line)){

    process_input_line( line, OUT) ;

/*
  Skip exclamation points in column 1 and blank lines
*/
    if ( OUT[0].empty() ){ continue ;}

    if ( OUT[0] == "hamilt" ){
/*
  Hamiltonian 
      Molecular/mol - 1
        no additional options needed
      Semi-Empirical/se - 2
        Specify the set of approximations
        Specify the parameter set
      Hubbard - 3
        Specify the interaction strength
        Specify the dimensions
      Pairing - 4
        Specify the interaction strength
        Specify the number of levels
*/

      if ( OUT[1] == "mol") {
        /* Molecular Electronic Hamiltonian */
        com.hamil(1) ;
      } else if ( OUT[1] == "se") {
        /* Semi-Empirical Hamiltonian */
        std::cout << "  NOT IMPLEMENTED " << std::endl ;
        com.hamil(1) ;
      } else if ( OUT[1] == "hub" || OUT[1] == "phub") {
        /* Hubbard */
        i_junk = 3 ;
        com.hamil(i_junk) ;
	if ( OUT[1].substr(0,1) == "p"){
	  /* Periodic */
          i_junk = 10 ;
        } else {
          i_junk = 0 ;
          }
        getline( jobfile, line) ;
        while( line.substr(0,4) != " end"){
          t_junk = 0 ;
  	  if ( line.substr(0,3) == " U "){
            d_junk = stod(line.substr(4)) ;
            com.setU( d_junk) ;
          } else if ( line.substr(0,3) == " nx"){
            i_junk = stoi(line.substr(4)) ;
            t_junk += i_junk ;
            com.hub_n( i_junk, 0) ;
          } else if ( line.substr(0,3) == " ny"){
            i_junk = stoi(line.substr(4)) ;
            t_junk += i_junk ;
            com.hub_n( i_junk, 1) ;
          } else if ( line.substr(0,3) == " nz"){
            i_junk = stoi(line.substr(4)) ;
            t_junk += i_junk ;
            com.hub_n( i_junk, 2) ;
          } else {
            qtzcntrl::shutdown( "Read Input: Unrecognized option in Hubbard") ;
            }
          getline( jobfile, line) ;
          }
        com.nbas( t_junk) ;
      } else if ( OUT[1] == "pair"){
        com.hamil( 4) ;
          getline( jobfile, line) ;
//          line.substr(0,3) == " U " ;
          i_junk = stoi(line.substr(7)) ;
          com.nalp( i_junk/2) ;
          com.nbet( i_junk/2) ;
          getline( jobfile, line) ;
//          line.substr(0,3) == " U " ;
          d_junk = stod(line.substr(5)) ;
          com.setU( d_junk) ;
          getline( jobfile, line) ;
//          line.substr(0,3) == " n " ;
          i_junk = stoi(line.substr(5)) ;
          com.nbas( i_junk) ;
          
      } else {
        com.hamil(1) ;
        }
/*

*/

    } else if ( OUT[0] == "method" ){
/*
  List of options :

    Single Reference Wavefunction
      This should be a string containing the possible symmetries to break
      and takes the form ABC where
        A - whether the wavefuntion is real or complex
          real/r
          complex/c

        B - how the spin symmetry should be treated
          restricted/r
          unrestricted/u
          genralized/g

        C - the type of wavefunction
          Hartree-Fock/HF
          Hartree-Fock-Bogliubov/HFB

    Projection methods
      P{..}|...>
        The options in the brakets {} are the types of projections
          N - Number projection
          T - Time Reversal
          S - Electron Spin
          K - Complex Conjugation

        The options in the ket will contain the type of reference 
        determinant to use listed below.
*/

/*
  Check for projection
*/
      if ( OUT[1].substr(0) == "p" ){

/*
  Get the projection options
*/
        first = OUT[1].find("{") ;
        last = OUT[1].find("}") ;
        s_junk = OUT[1].substr( first+1, last-2) ;
        boost::split( OUTx, s_junk, boost::is_any_of(",")) ;
        for ( it = OUTx.begin(); it != OUTx.end(); it++){
          strip_lower( *it) ;
          }

        for ( it = OUTx.begin(); it != OUTx.end(); it++){
          i_junk = 0 ;
          switch ( prjcase[*it]) {
            case 1 : /* number projection */
              i_junk |= 1 ;
              break ;
            case 2 : /* spin projection */
              i_junk |= 2 ;
              break ;
            case 3 : /* time reversal projection */
              i_junk |= 4 ;
              break ;
            case 4 : /* complex conjugation projection */
              i_junk |= 8 ;
              break ;
            default :
              qtzcntrl::shutdown( "Unrecognized Projection Symmetry") ;
              break ;
            } // End switch
          }

        com.methd( 100) ;
/*
  Parse the remaining projection parameters
*/
        first = OUT[1].find("|") ;
        last = OUT[1].find(">") ;
        OUT[1] = OUT[1].substr( first+1, last-2) ;
        } // End projection parsing

      if ( OUT[1].substr( 0, 2) == "rr" ) { com.methd(1) ;}
      else if ( OUT[1].substr(0,2) == "cr" ) { com.methd(2) ;}
      else if ( OUT[1].substr(0,2) == "ru" ) { com.methd(3) ;}
      else if ( OUT[1].substr(0,2) == "cu" ) { com.methd(4) ;}
      else if ( OUT[1].substr(0,2) == "rg" ) { com.methd(5) ;}
      else if ( OUT[1].substr(0,2) == "cg" ) { com.methd(6) ;} // End Symmetry parsing

      if ( OUT[1].substr(2) == "hf" ) {
        com.methd(10) ;
      } else if ( OUT[1].substr(2) == "hfb" ) {
        com.methd(20) ;
      } else {
        qtzcntrl::shutdown(" Unrecognized method ") ;
	}

      if ( com.methd()/10 != 1 ){
        getline( jobfile, line) ;

        process_input_line( line, OUTx) ;

        while( OUTx[0] != "end" ) {
          if ( OUTx[0] == "ngrid" ) {
            // Size of the number projection grid
            i_junk = stoi( OUTx[1]) ;
            com.ngrid(i_junk) ;
          } else if ( OUTx[0] == "mu" ) {
            // Initial guess for the chemical potential
            d_junk = stod( OUTx[1]) ;
            com.mu( d_junk) ;
          } else {
            std::cout << "No other options implemented currently" << std::endl ;
            }
          getline( jobfile, line) ;
          process_input_line( line, OUTx) ;
          } // End additional parameter parsing
        }
    } else if ( OUT[0] == "basis" ){
      com.bnam( OUT[1]) ;
    } else if ( OUT[0] == "guess" ){
      /*
        We have some limited options right now.  Its either true or false
      */
       com.ini_guess( true) ;
    }  else if ( OUT[0] == "geom" ) {
      t_c.clear() ;
      getline( jobfile, line) ;
      while( line.substr(0,5) != " end" ){
        at_count++ ;
        pos = line.find(delim) ;
        d_junk = stod(line.substr( 0, pos)) ;
        atnum.push_back( d_junk) ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[0] = stod(line.substr( 0, pos)) ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[1] = stod(line.substr( 0, pos)) ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[2] = stod(line.substr( 0, pos)) ;
        line.erase( 0, pos + delim.length()) ;
        getline( jobfile, line) ;
        t_c.push_back( tmp) ;
        }
      com.natm( at_count) ;
      com.setA( atnum) ;
      com.setC( t_c) ;
    } else if ( OUT[0] == "nspin" ) {
      t_c.clear() ;
      getline( jobfile, line) ;
      while( line.substr(0,5) != " end" ){
        pos = line.find(delim) ;
        d_junk = stod(line.substr( 0, pos)) ;
        /* 
          This holds the angular momentum of the atom.
        */
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[0] = stod(line.substr( 0, pos)) ;
        d_calc = tmp[0]*tmp[0] ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[1] = stod(line.substr( 0, pos)) ;
        d_calc += tmp[1]*tmp[1] ;
        line.erase( 0, pos + delim.length()) ;
        pos = line.find(delim) ;
        tmp[2] = stod(line.substr( 0, pos)) ;
        d_calc += tmp[2]*tmp[2] ;
        line.erase( 0, pos + delim.length()) ;
        /* 
          Normalize
        */
        if ( d_calc > 1.0e-8){
          d_junk /= std::sqrt( d_calc) ;
          tmp[0] *= d_junk ;
          tmp[1] *= d_junk ;
          tmp[2] *= d_junk ;
          }

        getline( jobfile, line) ;
        t_c.push_back( tmp) ;
        }
      com.setNS( t_c) ;
    } else if ( line.substr(0,7) == "options" ) {
/*
  Parse random options we may want to use
*/
      getline( jobfile, line) ;
      while( line.substr(0,5) != " end" ){
        if ( line.substr(0,7) == "mxscfit" ) {
          i_junk = stoi(line.substr(9)) ;
          com.mxscfit(i_junk) ;
        } else if ( line.substr(0,7) == "mxpnit" ) {
          i_junk = stoi(line.substr(9)) ;
          com.mxpnit(i_junk) ;
        } else if ( line.substr(0,7) == "print  " ) {
          i_junk = stoi(line.substr(9)) ;
          com.prt(i_junk) ;
        } else if ( line.substr(0,7) == "lshift " ) {
          d_junk = stod(line.substr(9)) ;
          com.lvlshft(d_junk) ;
          }
        getline( jobfile, line) ;
        }
    } else if ( line.substr(0,7) == "nalpha " ) {
      i_junk = stoi(line.substr(9)) ;
      com.nalp(i_junk) ;
    } else if ( line.substr(0,7) == "nbeta  " ) {
      i_junk = stoi(line.substr(9)) ;
      com.nbet(i_junk) ;
    }
  }

  jobfile.close() ;

  com.nele(com.nalp() + com.nbet()) ;
  read_input_time.end() ;

  return ;

}  ;

//void read_input( common& com, const std::string& inpfile){
// /* Read the input file and save the information into common. */
//  int i_junk, t_junk ;
//  int at_count = 0 ;
//  std::vector<double> tmp(3) ;
//  double d_junk ;
//  size_t pos = 0 ;
//  std::string line ;
//  std::string s_junk, s_j1, s_j2 ;
//  std::string delim = "," ;
//  std::string::size_type first, last ;
//  std::ifstream jobfile ;
//  std::vector<double> atnum ;
//  std::vector<std::vector<double>> t_c ;
//  std::unordered_map< std::string, int> prjcase ( { {"n", 1}, {"s", 2}, {"t", 3}, {"k", 4}} ) ;
//  time_dbg read_input_time = time_dbg("read_input") ;
//
//  jobfile.open( inpfile, std::ifstream::in ) ;
//
//  while( getline( jobfile, line)){
///*
//  Skip exclamation points in column 1
//*/
//    if ( line.substr(0,1) == "!" ){ continue ;}
//
//    if ( line.substr(0,7) == "hamilt " ){
//      s_junk = line.substr(9) ;
//      strip_lower( s_junk) ;
///*
//      Molecular/mol - 1
//      Semi-Empirical/se - 2
//      Hubbard - 3
//      Pairing - 4
//*/
//      if ( s_junk == "mol"){
//        /* Molecular Electronic Hamiltonian */
//        com.hamil(1) ;
//      } else if ( s_junk == "se"){
//        /* Semi-Empirical Hamiltonian */
//        std::cout << "  NOT IMPLEMENTED " << std::endl ;
//        com.hamil(1) ;
//      } else if ( s_junk.substr(0,3) == "hub" || s_junk.substr(0,4) == "phub"){
//        /* Hubbard */
//        i_junk = 3 ;
//        com.hamil(i_junk) ;
//	if ( s_junk.substr(0,1) == "p"){
//	  /* Periodic */
//          i_junk = 10 ;
//        } else {
//          i_junk = 0 ;
//          }
//        getline( jobfile, line) ;
//        while( line.substr(0,4) != " end"){
//          t_junk = 0 ;
//  	  if ( line.substr(0,3) == " U "){
//            d_junk = stod(line.substr(4)) ;
//            com.setU( d_junk) ;
//          } else if ( line.substr(0,3) == " nx"){
//            i_junk = stoi(line.substr(4)) ;
//            t_junk += i_junk ;
//            com.hub_n( i_junk, 0) ;
//          } else if ( line.substr(0,3) == " ny"){
//            i_junk = stoi(line.substr(4)) ;
//            t_junk += i_junk ;
//            com.hub_n( i_junk, 1) ;
//          } else if ( line.substr(0,3) == " nz"){
//            i_junk = stoi(line.substr(4)) ;
//            t_junk += i_junk ;
//            com.hub_n( i_junk, 2) ;
//          } else {
//            qtzcntrl::shutdown( "Read Input: Unrecognized option in Hubbard") ;
//            }
//          getline( jobfile, line) ;
//          }
//        com.nbas( t_junk) ;
//      } else if ( s_junk.substr(0,4) == "pair"){
//        com.hamil( 4) ;
//          getline( jobfile, line) ;
////          line.substr(0,3) == " U " ;
//          i_junk = stoi(line.substr(7)) ;
//          com.nalp( i_junk/2) ;
//          com.nbet( i_junk/2) ;
//          getline( jobfile, line) ;
////          line.substr(0,3) == " U " ;
//          d_junk = stod(line.substr(5)) ;
//          com.setU( d_junk) ;
//          getline( jobfile, line) ;
////          line.substr(0,3) == " n " ;
//          i_junk = stoi(line.substr(5)) ;
//          com.nbas( i_junk) ;
//          
//      } else {
//        com.hamil(1) ;
//        }
//    } else if ( line.substr( 0, 7) == "method " ){
///*
//  List of options :
//
//    Projection methods
//      P{..}|...>
//        The options in the brakets {} are the types of projections
//          N - Number projection
//          T - Time Reversal
//          S - Electron Spin
//          K - Complex Conjugation
//
//        The options in the ket will contain the type of reference 
//        determinant to use listed below.
//
//    Single Reference Wavefunction
//      This should be a string containing the possible symmetries to break
//      and takes the form ABC where
//        A - whether the wavefuntion is real or complex
//          real/r
//          complex/c
//
//        B - how the spin symmetry should be treated
//          restricted/r
//          unrestricted/u
//          genralized/g
//
//        C - the type of wavefunction
//          Hartree-Fock/HF
//          Hartree-Fock-Bogliubov/HFB
//*/
//      s_junk = line.substr(9) ;
//      strip_lower( s_junk) ;
///*
//  Check for projection
//*/
//      if ( s_junk.substr(0, 1) == "p" ){
///*
//  Get the projection options
//    The various combinations are stored as bitwise
//
//    N  - 0001 -> 1
//    S  - 0010 -> 2
//    K  - 0100 -> 4
//    T  - 1000 -> 8
//*/
//        /* Add 20 here to generate an initial guess */
//        com.methd(20) ;
//        s_junk.pop_back() ;
//        first = s_junk.find("{") ;
//        last = s_junk.find("}") ;
//        s_j1 = s_junk.substr( first+1, last-2) ;
//        s_junk.erase( 0, last+2) ;
//        do {
//          pos = s_j1.find(delim) ;
//          s_j2 = s_j1.substr(0, pos) ;
//          s_j1.erase( 0, pos + delim.length()) ;
//          i_junk = 0 ;
//          switch ( prjcase[s_j2]) {
//            case 1 : /* number projection */
//              i_junk |= 1 ;
//              break ;
//            case 2 : /* spin projection */
//              i_junk |= 2 ;
//              break ;
//            case 3 : /* time reversal projection */
//              i_junk |= 4 ;
//              break ;
//            case 4 : /* complex conjugation projection */
//              i_junk |= 8 ;
//              break ;
//            default :
//              qtzcntrl::shutdown( "Unrecognized Projection Symmetry") ;
//            } // End switch
//          } while ( pos != std::string::npos) ; // End While
//        
//          com.methd( i_junk*100) ;
//        getline( jobfile, line) ;
//        while( line.substr(0,5) != " end" ){
//          if ( line.substr( 0, 7) == " ngrid " ) { 
//            // Size of the number projection grid
//            i_junk = stoi(line.substr(9)) ;
//            com.ngrid(i_junk) ;
//            }
//          getline( jobfile, line) ;
//          } // End Project parameter parsing
//
//        } // End projection parsing
//
//      if ( s_junk.substr( 0, 2) == "rr" ) { com.methd(1) ;}
//      else if ( s_junk.substr(0,2) == "cr" ) { com.methd(2) ;}
//      else if ( s_junk.substr(0,2) == "ru" ) { com.methd(3) ;}
//      else if ( s_junk.substr(0,2) == "cu" ) { com.methd(4) ;}
//      else if ( s_junk.substr(0,2) == "rg" ) { com.methd(5) ;}
//      else if ( s_junk.substr(0,2) == "cg" ) { com.methd(6) ;} // End Symmetry parsing
//
//      if ( s_junk.substr(2) == "hf" ) { com.methd(10) ;}
//      else if ( s_junk.substr( 2) == "hfb" ) {
//        getline( jobfile, line) ;
//        s_junk = line.substr( 0, 7) ;
//        strip_lower( s_junk) ;
//        if ( s_junk == "mu") {
//          com.methd(20) ;
//        } else {
//          qtzcntrl::shutdown( "Unrecognized chemical potential") ;
//	  }
//        d_junk = stod(line.substr(9)) ;
//        com.mu( d_junk) ;
//        } // End HF/HFB parsing
//
//    } else if ( line.substr(0,7) == "basis  " ){
//      s_junk = line.substr(9) ;
//      com.bnam( s_junk) ;
//    }  else if ( line.substr(0,7) == "geom   " ) {
//      getline( jobfile, line) ;
//      while( line.substr(0,5) != " end" ){
//        at_count++ ;
//        pos = line.find(delim) ;
//        d_junk = stod(line.substr( 0, pos)) ;
//        atnum.push_back( d_junk) ;
//        line.erase( 0, pos + delim.length()) ;
//        pos = line.find(delim) ;
//        tmp[0] = stod(line.substr( 0, pos)) ;
//        line.erase( 0, pos + delim.length()) ;
//        pos = line.find(delim) ;
//        tmp[1] = stod(line.substr( 0, pos)) ;
//        line.erase( 0, pos + delim.length()) ;
//        pos = line.find(delim) ;
//        tmp[2] = stod(line.substr( 0, pos)) ;
//        line.erase( 0, pos + delim.length()) ;
//        getline( jobfile, line) ;
//        t_c.push_back( tmp) ;
//        }
//      com.natm( at_count) ;
//      com.setA( atnum) ;
//      com.setC( t_c) ;
//    } else if ( line.substr(0,7) == "options" ) {
///*
//  Parse random options we may want to use
//*/
//      getline( jobfile, line) ;
//      while( line.substr(0,5) != " end" ){
//        if ( line.substr(0,7) == "mxscfit" ) {
//          i_junk = stoi(line.substr(9)) ;
//          com.mxscfit(i_junk) ;
//        } else if ( line.substr(0,7) == "mxpnit" ) {
//          i_junk = stoi(line.substr(9)) ;
//          com.mxpnit(i_junk) ;
//        } else if ( line.substr(0,7) == "print  " ) {
//          i_junk = stoi(line.substr(9)) ;
//          com.prt(i_junk) ;
//        } else if ( line.substr(0,7) == "lshift " ) {
//          d_junk = stod(line.substr(9)) ;
//          com.lvlshft(d_junk) ;
//          }
//        getline( jobfile, line) ;
//        }
//    } else if ( line.substr(0,7) == "nalpha " ) {
//      i_junk = stoi(line.substr(9)) ;
//      com.nalp(i_junk) ;
//    } else if ( line.substr(0,7) == "nbeta  " ) {
//      i_junk = stoi(line.substr(9)) ;
//      com.nbet(i_junk) ;
//    }
//  }
//
//  jobfile.close() ;
//
//  com.nele(com.nalp() + com.nbet()) ;
//  read_input_time.end() ;
//
//  return ;
//
//}  ;

bool open_text( std::ofstream& F_OUT, int cntl, const std::string& filename) {
/*
  cntl :
    0 - (default) write to a new file or overwrite one that already exists
    1 - Check if a file exists and append.  If it does not exist create it.
*/
  const std::string wfnIO = filename + ".txt" ;
  bool ioerr = false ;
  struct stat buf ;

  if ( cntl == 0 ){
    try {
      F_OUT.open( wfnIO, std::ofstream::out) ;
    } catch ( std::fstream::failure& e) {
      std::cout << e.what() << std::endl ;
      ioerr = true ;
    } catch (...) {
      std::cout << " Default excepetion inside write_to_bin " << std::endl ;
      ioerr = true ;
      }
  } else if ( cntl == 1 ){
    /* Check if the file exists */
    if ( stat(wfnIO.c_str(), &buf) == 0) {
      try {
        F_OUT.open( wfnIO, std::ofstream::app | std::ofstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
    } else {
      try {
        F_OUT.open( wfnIO, std::ofstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
    }
  } else {
    qtzcntrl::shutdown( "Unrecognized Option (open_binary)") ;
    }

  return ioerr ;

} ;

bool open_binary( std::ofstream& F_OUT, int cntl) {
/*
  cntl :
    0 - (default) write to a new file or overwrite one that already exists
    1 - Check if a file exists.  If not create it.
*/
  const std::string wfnIO = "qtz.wfn.bin" ;
  bool ioerr = false ;
  struct stat buf ;

  if ( cntl == 0 ){
    try {
      F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
    } catch ( std::fstream::failure& e) {
      std::cout << e.what() << std::endl ;
      ioerr = true ;
    } catch (...) {
      std::cout << " Default excepetion inside write_to_bin " << std::endl ;
      ioerr = true ;
      }
  } else if ( cntl == 1 ){
    /* Check if the file exists */
    if ( stat(wfnIO.c_str(), &buf) == 0) {
      try {
        F_OUT.open( wfnIO, std::ofstream::app | std::ofstream::binary | std::ofstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
    } else {
      try {
        F_OUT.open( wfnIO, std::ofstream::binary | std::ofstream::out) ;
      } catch ( std::fstream::failure& e) {
        std::cout << e.what() << std::endl ;
        ioerr = true ;
      } catch (...) {
        std::cout << " Default excepetion inside write_to_bin " << std::endl ;
        ioerr = true ;
        }
    }
  } else {
    qtzcntrl::shutdown( "Unrecognized Option (open_binary)") ;
    }

  return ioerr ;

} ;

bool open_binary( std::ifstream& F_IN) {
/*
  (default) Check if a file exists and if it does open it
*/
  const std::string wfnIO = "qtz.wfn.bin" ;
  bool ioerr = false ;
  struct stat buf ;

  /* Check if the file exists */
  if ( stat(wfnIO.c_str(), &buf) == 0) {
    try {
      F_IN.open( wfnIO, std::ofstream::binary | std::ofstream::in) ;
    } catch ( std::fstream::failure& e) {
      std::cout << e.what() << std::endl ;
      ioerr = true ;
    } catch (...) {
      std::cout << " Default exception inside write_to_bin " << std::endl ;
      ioerr = true ;
      }
  } else {
    qtzcntrl::shutdown( "FIle does not exist (open_binary)") ;
    }

  return ioerr ;

} ;

void process_input_line( std::string& raw, std::vector<std::string>& OUT) {
  std::vector<std::string>::iterator it ;
/*
  Remove anything after an exclamation mark
*/

  boost::split( OUT, raw, boost::is_any_of("!")) ;
  raw = OUT[0] ;

  boost::split( OUT, raw, boost::is_any_of(":")) ;
  for ( it = OUT.begin(); it != OUT.end(); it++){
    strip_lower( *it) ;
    }

  } ;

void strip_lower( std::string& s) {
/*
  Wrap together whitespace trimming and lower casing for input file parsing
*/
  std::string::iterator end_pos ;
  end_pos = std::remove( s.begin(), s.end(), ' ') ;
  s.erase( end_pos, s.end()) ;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower) ;

  return ;

} ;

template <class matrix> 
void write_eigen_bin (const matrix& m, std::ofstream& F_OUT) {
  typename matrix::Index rows=m.rows(), cols=m.cols();

  F_OUT.write( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;

  return ;

}

template void write_eigen_bin(const Eigen::VectorXd&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::VectorXcd&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::MatrixXd&, std::ofstream&) ;
template void write_eigen_bin(const Eigen::MatrixXcd&, std::ofstream&) ;

template <class matrix> 
void read_eigen_bin (const matrix& m, std::ifstream& F_IN) {
  typename matrix::Index rows=m.rows(), cols=m.cols();

  F_IN.read( (char*) m.data(), rows*cols*sizeof(typename matrix::Scalar)) ;

  return ;

}

template void read_eigen_bin(const Eigen::VectorXd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::VectorXcd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXd&, std::ifstream&) ;
template void read_eigen_bin(const Eigen::MatrixXcd&, std::ifstream&) ;

template <class matrix>
void print_mat( const matrix& o, std::string h){
  int i = 0 ;
  typename matrix::Index cols = o.cols(), rows = o.rows() ;
  Eigen::IOFormat matprt(5, 0, "  ", "\n", "| ", "|", "") ;

  if ( ! h.empty()){
    std::cout << h << std::endl ;
    std::cout << "---- ---- ----" << std::endl ;
    }

  if ( cols <= 5){
      std::cout << o.format( matprt) << std::endl ;
  } else {
    do {
      std::cout << o.block( 0, i, rows, 5).format( matprt) << std::endl << std::endl ;
      i += 5 ;
      cols -= 5 ;
      } while ( cols > 5) ;
      std::cout << o.block( 0, i, rows, cols).format( matprt) << std::endl ;
    }

  return ;

}

template void print_mat(const Eigen::VectorXd& , std::string) ;
template void print_mat(const Eigen::VectorXcd& , std::string) ;
template void print_mat(const Eigen::MatrixXd& , std::string) ;
template void print_mat(const Eigen::MatrixXcd& , std::string) ;
template void print_mat(const Eigen::Ref<Eigen::VectorXd>& , std::string) ;
template void print_mat(const Eigen::Ref<Eigen::VectorXcd>& , std::string) ;
template void print_mat(const Eigen::Ref<Eigen::MatrixXd>& , std::string) ;
template void print_mat(const Eigen::Ref<Eigen::MatrixXcd>& , std::string) ;

