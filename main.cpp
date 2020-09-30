/*
    Quetzalcoatl - electronic structure package
    Copyright (C) 2020  Kyle Throssell
    email :: ktthross@gmail.com
*/
#include <iostream>
#include <libint2.hpp>
#include "qtzctl.h"
#include "qtzio.h"
#include <string>
#include <vector>

void quetzalcoatl(const bool print){
/* Title Card */
  if (print) {
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
    }
  return ;
}

int main(int argc, char *argv[]) {
  
  quetzalcoatl(false);
  QtzInput parser(argc, argv);
  parser.parse_input();
  QtzControl qtz_control = parser.control();
  std::vector<libint2::Atom> atoms = parser.atoms();
  if (qtz_control.directive == "wavefunction"){
    std::cout << "Solving for a wavefunction" << std::endl;
    }
  /*
   * Figure out what we need to do.  Let's start with solving for wavefunction
   * */

  return 0 ;

} /* End Quetzacoatl */

