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
#include "qtzio.h"
#include "sladet.h"
#include "solver.h"
#include "tei.h"
#include "time_dbg.h"
#include "util.h"

void quetzalcoatl( void){
/* Title Card */
  std::cout << std::endl ;
  std::cout << "________                 __                 .__                      __  .__   " << std::endl ;
  std::cout << "\\_____  \\  __ __   _____/  |______________  |  |   ____  _________ _/  |_|  |  " << std::endl ;
  std::cout << " /  / \\  \\|  |  \\_/ __ \\   __\\___   /\\__  \\ |  | _/ ___\\/  _ \\__  \\\\   __\\  |  " << std::endl ;
  std::cout << "/   \\_/.  \\  |  /\\  ___/|  |  /    /  / __ \\|  |_\\  \\__(  <_> ) __ \\|  | |  |__" << std::endl ;
  std::cout << "\\_____\\ \\_/____/  \\___  >__| /_____ \\(____  /____/\\___  >____(____  /__| |____/" << std::endl ;
  std::cout << "       \\__>           \\/           \\/     \\/          \\/          \\/           " << std::endl ;
  std::cout << std::endl ;
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
  std::cout << "            oy- ```:ooooooossoo++o+oso+oss++/-:---.......--`  .../oo+ss+++++++" << std::endl ;
  std::cout << "            o+ss+-`..`-/+osooooooooosys+-:+oo/:-.-:-----.```....:--/ooooo+++++" << std::endl ;
  std::cout << "           ohsy++osso++//:::---:///++ooss+ossooossoo/://::---:+ohdNmyosso+++++" << std::endl ;
  std::cout << "                ++++++++ooooo+++++++++++++++++++oosooosssssssoo+++++++++++++++" << std::endl ;
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
  std::vector<tei> trnint ;

  int job=0  ;
  int i, j, k, l ;
  int nbas, natm ;
  int iopt ;
  double cx = 0.0e0 ;
  double cy = 0.0e0 ;
  double cz = 0.0e0 ;
  double r = 0.0e0 ;
  double r2 = 0.0e0 ;
  double n_rep = 0.0e0 ;
  double val = 0.0e0 ;
  cd ejunk ;
  cd ojunk ;
  std::vector<std::string> wfn_vec ;
  std::string trden="tden.rwf" ;
  std::string fokmat="fmat.rwf" ;
  Eigen::MatrixXd c ;
  Eigen::VectorXd a ;
  Eigen::MatrixXd S ;
  Eigen::MatrixXd T ;
  Eigen::MatrixXd V ;
  Eigen::MatrixXcd cV ;
  basis_set b ;
  sladet< double, Eigen::Dynamic, Eigen::Dynamic> w ;
  std::ofstream tstfile ; 
  std::ifstream tstfe ; 
  time_dbg quetz_time = time_dbg("Quetzalcoatl") ;

  /* File reading and header variables. */
  std::stringstream ss ;
  std::string inpfile ;

  quetzalcoatl() ;

  ss << argv[1] ;
  ss >> inpfile ;
  read_input( com, inpfile) ;

  /* Step 1 :
     Build the relevant data in memory.
     SCF routines
     real/complex reastricted
                  unrestricted
                  generalized  */
  
  natm = com.natm() ;
  c.resize( natm, 3) ;
  a.resize( natm) ;
  c = com.getC() ;
  a = com.getA() ;
  for (int i=0; i < natm; i++){
    for (int j=i+1; j < natm; j++){
      cx = c( j, 0) - c( i, 0) ;
      cy = c( j, 1) - c( i, 1) ;
      cz = c( j, 2) - c( i, 2) ;
      r2 = pow( cx, 2.0) + pow( cy, 2.0) + pow( cz, 2.0) ;
      n_rep += a(i)*a(j)/sqrt(r2) ;
      }
    }
 
  com.nrep( n_rep ) ;

  b = build_basis( com.bnam(), a, c) ;
  nbas = b.nbas ;
  com.nbas( nbas) ;
  S.resize( nbas, nbas) ;
  T.resize( nbas, nbas) ;
  V.resize( nbas, nbas) ;
  cV.resize( nbas, nbas) ;
  ao_overlap( com.natm(), b, S) ;
  com.setS( S) ;
  // Find the orthogonalizing routine
  canort( S, cV, nbas) ;
  T = -cV.real() ;
  com.setXS( T) ;
  ao_kinetic( com.natm(), b, T) ;
  ao_eN_V( com.natm(), b, c, a, V) ;
  S = T + V ;
  com.setH( S) ;
  T.resize( 0, 0) ;
  V.resize( 0, 0) ;
  S.resize( 0, 0) ;
  list_ao_tei( com.natm(), b, intarr) ; 

  scf_drv( com, intarr, com.methd()) ;
  mo_integrals( com, intarr) ;

  intarr.clear() ;
  V.resize( 0, 0) ; 
  T.resize( 0, 0) ; 
  S.resize( 0, 0) ;  
  quetz_time.end() ;

  return 0 ;

} /* End Quetzacoatl */

