/** 
 * @file temp.cc
 * Calculate the temperature for different sets of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor temp
 * @section temp Calculates the temperature for different sets of atoms
 * @author nb
 * @date 1.6.2012
 *
 * Program temp will calculate the temperatur for different sets of atoms,
 * as specified by the atomspecifier(s). 
 * 
 * You can find the number of degree of freedoms for a temperature group in the
 * md++ output file under "DEGREES OF FREEDOM" -> "DOF"
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>]
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier Atoms sets &gt; </td></tr>
 * <tr><td> \@dofs</td><td>&lt;@ref degrees of freedom &gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;velocity trajectory files&gt; </td></tr>
 * </table>
 *
 * @sa AtomSpecifier
 *   
 * Example:
@verbatim
 temp
    @topo       ex.top
    @time       0 0.1
    @atoms      1:a
    @dofs       26
    @traj       ex.trv

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"

#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/Temperature.h"

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace std;
using namespace gmath;

int main(int argc, char **argv){
  
  Argument_List knowns; 
  knowns << "topo" << "traj" << "atoms" << "time" << "dofs";

  string usage = "# " + string(argv[0]);
  usage += "# You can find the number of degree of freedoms for a temperature group in the\n";
  usage += "# md++ output file under \"DEGREES OF FREEDOM\" -> \"DOF\".\n";
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@atoms      <atoms to consider for the temperature calculations>\n";
  usage += "\t@dofs       <the number of degree of freedoms for these atoms>\n";
  usage += "\t@traj       <trajectory files>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // Check the number of arguments
    bool error = false;
    if (!(args.count("atoms") > 0)){
      cerr << "No atoms specified with '@atoms'!\n";
      error = true;
    }
    if (!(args.count("dofs") > 0)){
      cerr << "The number of degree of freedoms is not specified with '@dofs'!\n";
      error = true;
    }
    if (!(args.count("traj") > 0)){
      cerr << "No trajectory files specified with '@traj'!\n";
      error = true;
    }
    if (!(args.count("topo") > 0)){
      cerr << "No No topology specified with '@topo'!\n";
      error = true;
    }
    if (args.count("atoms") != args.count("dofs")){
      cerr << "There are not as many degree of freedom numbers as atom specifiers!\n";
      error = true;
    }
    if (error){
      cerr << "\nErrors during arguments checking!\n";
      exit(1);
    }
  
    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // Temperature
    bool has_solvent = false;
    cout << "# Time ";
    int temp_group = 1;
    vector<Temperature> temps;
    
    Arguments::iterator as_i = args.lower_bound("atoms"),
            as_to = args.upper_bound("atoms");
    Arguments::iterator dofs_i = args.lower_bound("dofs");
    
    for (; as_i != as_to; as_i++, dofs_i++){
      stringstream ss_dof;
      ss_dof << dofs_i->second;
      double dofs;
      ss_dof >> dofs;
      temps.push_back(Temperature(AtomSpecifier(sys, as_i->second), dofs));
      cout << "<Temperature Group " << temp_group++ << ">";
      if (as_i->second.find('s') < as_i->second.size()){
        has_solvent = true;
      }
    }
    cout << endl;
    
    if (has_solvent){
      cout << "# Solvent specified!" << endl;
      cout << "# It might be better to turn the solvents into solutes in the topology!" << endl;
    }
    
    
    // Process the files and generate the temperature
    Arguments::iterator file_i = args.lower_bound("traj"),
            file_to = args.upper_bound("traj");
    for (; file_i != file_to; file_i++){
      InG96 ic;
      ic.open(file_i->second);
      while(!ic.eof()){
        ic.select("ALL");
        ic >> sys >> time;
        if (!sys.hasVel){
          cerr << "No Velocity read!\n";
          exit(1);
        }

        cout << time;
        vector<Temperature>::iterator temp_i = temps.begin(),
            temp_to = temps.end();
        for (; temp_i != temp_to; temp_i++){
          cout << " \t" << temp_i->temperature(sys);
        }
        cout << endl;
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
