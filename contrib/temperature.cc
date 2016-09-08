#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutPdb.h"

#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "delaunay.cc"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;


//JAG 2.6.2013
/* Temperature calculation cnf,
*/ 

const double KB=0.008314462;
int main(int argc, char **argv)  {
  
  
  Argument_List knowns;
  knowns << "topo" << "atoms"<< "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@traj   <trajectory files>\n";
    

  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    
    //read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    InG96 ic;
    if (args.count("traj") > 0) {
      ic.open(args.lower_bound("traj")->second);
      ic.select("ALL");
      ic >> refSys;
      ic.close();
    } else {
      throw gromos::Exception("frameout", "no trajectory specified (@traj)");
    }

   
    // The current system
    System sys(refSys);
    
    // get temperature atoms
    string spec;
    AtomSpecifier atoms_temp(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
	//cout<<"hola"<<endl;
        spec = iter->second.c_str();
	//cout<<spec<<endl;
        atoms_temp.addSpecifier(spec);
      }
    }
    if ( !atoms_temp.size() )
      throw gromos::Exception("temperature",
			      "No atoms to calculate temperature!");
    
    double kinetic=0.0;
    double temperature=0.0;
    // loop over all trajectories 
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to   = args.upper_bound("traj");
    for (; iter != to; iter++) {
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()) {
	ic >> sys >> time;
	for(int i = 0; i < atoms_temp.size(); i++) {
	  kinetic+=0.5*atoms_temp.mass(i)*(atoms_temp.vel(i)[0]*atoms_temp.vel(i)[0]);
	  kinetic+=0.5*atoms_temp.mass(i)*(atoms_temp.vel(i)[1]*atoms_temp.vel(i)[1]);
	  kinetic+=0.5*atoms_temp.mass(i)*(atoms_temp.vel(i)[2]*atoms_temp.vel(i)[2]);
	}
	kinetic/=1000;
	temperature=(2*kinetic)/(3*atoms_temp.size()*KB);
      }
      
    }
    
    cout<<"Temperature of "<<spec<<" = "<<temperature<<endl;
    
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
