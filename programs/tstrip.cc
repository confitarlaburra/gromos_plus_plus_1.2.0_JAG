/**
 * @file tstrip.cc
 * Removes solvent coordinates from a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor tstrip
 * @section tstrip Removes solvent coordinates from a trajectory
 * @author @ref rb
 * @date 11-6-07
 *
 * Program tstrip removes all solvent coordinates from a (list of) trajectory
 * file(s) for ease of later analysis. Note that program @ref filter 
 * (see section V-4.2) captures the functionality of tstrip as well.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input trajectory file(s)&gt; </td></tr>
 * <tr><td> [\@nthframe</td><td>&lt;write every nth frame (default: 1)&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  tstrip
    @topo      ex.top
    @traj      ex.tr
    @nthframe  3
    [@time     0 2]
 @endverbatim
 *
 * <hr>
 */

#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/groTime.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;



int main(int argc, char **argv){
  Argument_List knowns;
  knowns << "topo" << "traj" << "nthframe" << "time";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@traj <input trajectory files>\n";
  usage += "\t[@nthframe <write every nth frame> (default: 1)]\n";
  usage += "\t[@time      <time and dt>]\n";

  try{
    Arguments args(argc, argv, knowns, usage);
    
    utils::Time time(args);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    InG96 ic;
    OutCoordinates *oc;
    oc = new OutG96();

    int nthFrame = args.getValue<int>("nthframe", false, 1);

    // loop over all trajectories
    bool isFirstTraj = true;
    int skipFrame = 0;
    oc->open(cout);     
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);

      if (isFirstTraj){
        oc->writeTitle(ic.title());
        isFirstTraj = false;
      }

      // loop over all frames
      while(!ic.eof()){
        ic >> sys >> time;
        if (! skipFrame) {
          oc->writeTimestep(time.steps(), time.time());
          *oc << sys;
        }
        skipFrame++;
        skipFrame %= nthFrame;
      }    
    
      ic.close();
    }
    oc->close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
