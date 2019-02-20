/* Primary routine for Quetzalcoatl */
#include "basis.h"
#include "constants.h"
#include <complex>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <sstream>
#include "binio.h"
#include "common.h"
#include "evalm.h"
#include "hfrout.h"
#include "integr.h"
#include "obarasaika.h"
#include "postscf.h"
#include "project.h"
#include "nbodyint.h"
#include "qtzio.h"
#include "qtzcntrl.h"
#include "solver.h"
#include "tei.h"
#include "time_dbg.h"
#include "util.h"
#include "r12.h"
#include "wfn.h"
#include "staging.h"

/*
    Quetzalcoatl - electronic structure package
    Copyright (C) 2018  Kyle Throssell

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    email :: ktthross@gmail.com
*/

void quetzalcoatl( int im){
/* Title Card */
  std::cout << std::endl ;
  std::cout << "________                 __                 .__                      __  .__   " << std::endl ;
  std::cout << "\\_____  \\  __ __   _____/  |______________  |  |   ____  _________ _/  |_|  |  " << std::endl ;
  std::cout << " /  / \\  \\|  |  \\_/ __ \\   __\\___   /\\__  \\ |  | _/ ___\\/  _ \\__  \\\\   __\\  |  " << std::endl ;
  std::cout << "/   \\_/.  \\  |  /\\  ___/|  |  /    /  / __ \\|  |_\\  \\__(  <_> ) __ \\|  | |  |__" << std::endl ;
  std::cout << "\\_____\\ \\_/____/  \\___  >__| /_____ \\(____  /____/\\___  >____(____  /__| |____/" << std::endl ;
  std::cout << "       \\__>           \\/           \\/     \\/          \\/          \\/           " << std::endl ;
  std::cout << std::endl ;
  std::cout << " Quetzalcoatl  Copyright (C) 2018  Kyle Throssell " << std::endl ;
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY " << std::endl ;
  std::cout << std::endl ;
  if ( im == 1){
    std::cout << "        ooo+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl ;
    std::cout << "       o/`s++++++++++++++++++++++++++o++++++++++++oo+++//++o++++++++++++++" << std::endl ;
    std::cout << "      s:  o++++++++++++++++++++++++++oo+++++++++++o//.`    `/yo+++++++++++" << std::endl ;
    std::cout << "     s/   +++++++++++++++++++++++++++o/++++++++++++smhdy- .  `oyo+++++++++" << std::endl ;
    std::cout << "    oo /  ++++++++++++++so++++++++++++o./o++++++++++odymdyy`   /do++++++++" << std::endl ;
    std::cout << "    y`/:  o+++++++++++++oo+/o+++++++++oy..+o+++++++++odhdNmm/   -ho+++++++" << std::endl ;
    std::cout << "    d`s`  s++++++++++++++oso.-/oo++++++oo-`:o+++++++++ohmNNNy+   -do++++++" << std::endl ;
    std::cout << "   oss/:  s++++++++++++++++oy+``-++oo+++os/../o++++++++osmNmNm+   /mo+++++" << std::endl ;
    std::cout << "   s-y+.  s+++++++++++++++++oysy   `.://+o++-.-+o++++++++oyhMMd/   dh+++++" << std::endl ;
    std::cout << "   y`sy-  s+++++++++++++++++++oso---: -:-.+sso-.:+o++++++++ohhNh:  -No++++" << std::endl ;
    std::cout << "   s +s-  o+++++++++++++++++++++ooooooo+//y`/:---/yoo++++++++oyyh   dd++++" << std::endl ;
    std::cout << "   s ++`  ys+++++++++++++++++++++osy/:ooooyoo+:/:---/o++++++++oym.  +ds+++" << std::endl ;
    std::cout << "   s +o`  hh+++++++++++++++++oss+:..+ossoo/+:.-\\\\:./..o+++++++++m.  .oyo++" << std::endl ;
    std::cout << "   o`+o/  hNs++++++++++++++os+.`    .:+o:-/'\\{/0}}:/+.+o+++++++oy   :o:h++" << std::endl ;
    std::cout << "   o-//s  sMho+++++++++++os+.         ` ----s/:/ `\\```/o+++++++o+   sy:yo+" << std::endl ;
    std::cout << "    o.:y: /Myy+++++++++oo/-        `.-/+-://yhso:oo/-.-:o++++++y: .:s+/:h+" << std::endl ;
    std::cout << "    o/:h: -Ndy+++++++os:-.        /yosyoooooo+osyyyss++.++++++os`.ohy./ do" << std::endl ;
    std::cout << "     yos.`-ymms+++++oo` /       `ohyoo++++++++++ooyyooooo++++oy- /s+/`- ss" << std::endl ;
    std::cout << "     sds/ .mymdo+++oy` .-     `oyso+++++ooosooo+++++++++++++oy- .+`- -/`/d" << std::endl ;
    std::cout << "     oyyo  NymNy+++s:  /`     yho++++ooo+-.-.-:/+oooo++++++os- /+: :`:+-:m" << std::endl ;
    std::cout << "      sso` ymMNhs++h   /      dyosssso:`  .:--. `..-/+oooos:`  y.  +-.s-sd" << std::endl ;
    std::cout << "      s:++.`yMMhhoos-::..  `  s/.-oyd-  .:o/:/--///:````..`  ./s/ .o`-y/my" << std::endl ;
    std::cout << "      oo-s/-`/dMsy+```-..  +  `   /o/   :oso-./-.:`+`+.:`..:/o/`s.oo odydy" << std::endl ;
    std::cout << "       so:///``://  ./-+-`.+      -sy- `s++-/+:y.h-y-h /.--:`.s h-h`/+d-ms" << std::endl ;
    std::cout << "       ys-+:ss-   `:-:+.s`-+-   -. o- --yd-oss:.:-.:.o-o//+/ +N:ho/-yo./d+" << std::endl ;
    std::cout << "       ooo/oss++.`+-.`//++..::` .--:. ho+d`m//. h. .`/.`:::+//oo++:d:. ss+" << std::endl ;
    std::cout << "        oooy:sos//s..:/+:-+:`..    .:.d-+::m-/s do  ++. / `s/ /-o/:s:-.y++" << std::endl ;
    std::cout << "         ooo+o-s///-:--yoy+/+`      ./s `/./yyN.d+- oss-+:o-/oh+- ooo`so++" << std::endl ;
    std::cout << "           oo++s/o../+o+---:o/-` ``` `+s:+yoNd++o-s-ym:`-y:/hmM. .ds./s+++" << std::endl ;
    std::cout << "             oooyyo+::/:-``--...-` `.` `:/+///`  oo-hh.odssdMNN  yo.:s++++" << std::endl ;
    std::cout << "                odmhs+:///:.`--:-/+  `.`     `:+-shysssosymmsyd`oo./o+++++" << std::endl ;
    std::cout << "                 oymmmh+oyhhsso///s.   ```     `:/ho++odNmyshdysy:oo++++++" << std::endl ;
    std::cout << "                 oooyddhyydNMMMNdyoy`    .        .hoomdyhNMMsmyoso+++++++" << std::endl ;
    std::cout << "               oo/-.s.     `:odmmNMNo`````.        oshdNMMNddhmyo+++++++++" << std::endl ;
    std::cout << "              oy/   +/`     ``-syyys/--:::/:-`     ++.+my+/dmho+++++++++++" << std::endl ;
    std::cout << "              s-/`   :oo:`-/:://..          `--/++/:  `.s-yho++++++o++++++" << std::endl ;
    std::cout << "             oo.-+.    `:+so++-``              `+y`  `` -do+++++++os++++++" << std::endl ;
    std::cout << "            ooy`  .-.       `-:--::/:.`           ````   yo+++++osy+++++++" << std::endl ;
    std::cout << "          oo+:s-  ..--.`            `.++/-`   `         `y+++o+osys+++++++" << std::endl ;
    std::cout << "        os+.- .hhs-`   `-:.`         `y``-/o++..---.-..:so+o+:--syo+++++++" << std::endl ;
    std::cout << "        y/ -  .Nhhdhy+::-  `.+...-:-:os/+ooo+..``.--:+///--.`./syooo++++++" << std::endl ;
    std::cout << "        oy- ```:ooooooossoo++o+oso+oss++/-:---.......--`  .../oo+ss+++++++" << std::endl ;
    std::cout << "        o+ss+-`..`-/+osooooooooosys+-:+oo/:-.-:-----.```....:--/ooooo+++++" << std::endl ;
    std::cout << "       ohsy++osso++//:::---:///++ooss+ossooossoo/://::---:+ohdNmyosso+++++" << std::endl ;
    std::cout << "            ++++++++ooooo+++++++++++++++++++oosooosssssssoo+++++++++++++++" << std::endl ;
  } else if ( im == 2){
  } else if ( im == 3){
    }
  return ;
}

int main(int argc, char *argv[]) {
/*
 * This is my own implementation of various electronic stucture methods.  
 * There is lots of work to do and improvements to be made and the purpose
 * is supposed to be minamly pedagoical as well as producing results. 
 * Currently everything is done in memory.
 * To Do :: 
 *   - Parse options and set route through program.
 *   -implement stability 
 *   -matrix elements between determinants
 *   -Eventually, I will need to write a routine which grabs a chunk of memory,
 *   rather than constantly allocating and dealllocating.  For now, I need things
 *   that work and produce useful results.  
 *   - Debug spin rotations on a determinant.
 *   - Interface with non-gaussian routines for matrix elements
 *
 * */

/* 
  job_file - file to read information about the job
  readfile - name of the job file being read in 
  line - holds information being read in
  job_dim - A class to carry relevant dimensioning information through
	the program.
  
  ovl - matrix containing the overlp matrix
  ham - matrix containing the core hamiltonian
  xmat - matrix containing the transformation to an orthogonal ao basis.
  */
  common com = common() ;
  std::vector<tei> intarr ;

  int quetz_image = 0 ;
  cd ejunk ;
  cd ojunk ;
  std::vector<std::string> wfn_vec ;
  wfn< double, Eigen::Dynamic, Eigen::Dynamic> w ;
  std::ofstream tstfile ; 
  std::ifstream tstfe ; 
  std::srand((unsigned int) time(0)) ;
  time_dbg quetz_time = time_dbg("Quetzalcoatl") ;

  /* File reading and header variables. */
  std::stringstream ss ;
  std::string inpfile ;

  ss << argv[1] ;
  ss >> inpfile ;
  if ( inpfile == "Q" ){ 
    quetz_image = 1 ;
    ss.clear() ;
    ss << argv[2] ;
    ss >> inpfile ;
  } else if ( inpfile == "w" ){
    quetz_image = 2 ;
  } else if ( inpfile == "c" ){
    quetz_image = 3 ;
    }


  quetzalcoatl( quetz_image) ;

  read_input( com, inpfile) ;

/*
  Set up the appropriate Hamiltonian
*/

  if ( com.hamil() == 1 || com.hamil() == 2 ){

/*
  Molecular Hamiltonian
*/
//    std::vector<tei> intarr ;
//    molecular_hamiltonian( com, intarr) ;
    molecular_hamiltonian( com) ;

  } else if ( com.hamil() == 3 ){

/*
  Hubbard
*/
     ;
  } else if ( com.hamil() == 4 ){

/*
  Pairing
*/
  pairing_hamiltonian( com) ;

  } else {
    qtzcntrl::shutdown( "Unrecognized Hamiltonian in Main." ) ;
    }

/*
  Follow the appropriate path through the program
*/

  if ( true){

/*
  We always do a mean-field calculation so let's not worry about logic here quite yet
*/
    scf_drv( com) ;
    prj_drv( com) ;

    }

/*
  int cghfxx = 11 ;
  scf_drv( com, intarr, cghfxx) ;
  cghfxx = 21 ;
  scf_drv( com, intarr, cghfxx) ;


  if ( com.methd() / 100 != 0 ){
    int rrphfb = 121 ;
    prj_drv( com, intarr, rrphfb) ;
    }
*/
  quetz_time.end() ;

  return 0 ;

} /* End Quetzacoatl */

