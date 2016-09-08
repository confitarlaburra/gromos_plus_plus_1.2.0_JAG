#ifdef OMP
#include <omp.h>
#endif

#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "../src/args/Arguments.h"
//#include "../src/args/BoundaryParser.h"
//#include "../src/args/GatherParser.h"
//#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
//#include "../src/fit/RotationalFit.h"
//#include "../src/fit/TranslationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/args/OutformatParser.h"
//#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
//#include "../src/gio/OutPdb.h"

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

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;


//JAG 19.11.2013
/* Gathering for nanotubes,
*/ 




int main(int argc, char **argv)  {
  
  
  Argument_List knowns;
  knowns << "topo"  <<"atoms" << "traj" << "center" << "time" << "outformat" << "threads"<<"first"<<"last";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@center  <index of center solute atom>\n";
  usage += "\t@atoms  <atoms to gather>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@outformat  <output coordinates format>\n";
  usage += "\t@threads  <number of threads for openmp>\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @first argument
    int first = args.getValue<int>("first");
    
    // get the @last argument
    int last = args.getValue<int>("last");
        
    // get the @time argument, create time object
    utils::Time time(args);
    
    // get simulation time either from the user or from the files
    bool usertime=false;
    
    if (args.count("time") > 0) {
      usertime=true;
    }

    int threads;
    if (args.count("threads") > 0) 
      threads=args.getValue<int>("threads");
    else
      threads=1;

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
      throw gromos::Exception("gatherNT", "no trajectory specified (@traj)");
    }

   
    // The current system
    System sys(refSys);
    
    // get gather  atoms
    string spec;
    AtomSpecifier atoms_gather(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        spec = iter->second.c_str();
        atoms_gather.addSpecifier(spec);
      }
    }
    if ( !atoms_gather.size() )
      throw gromos::Exception("gather",
			      "No atoms to be gathered!");
    
    int index_center=args.getValue<int>("center");
    Vec max(0.0,0.0,0.0);
    Vec min(0.0,0.0,0.0);
    Vec vec_box;

    
    // output
    string ext;
    OutCoordinates & oc = *OutformatParser::parse(args, ext);
    oc.open(cout);
    std::ostringstream title;

    title << "Gathered trajectory. Keeping" << endl;
    if (atoms_gather.size()) {
      vector<string> s = atoms_gather.toString();
      title << "* atoms ";
      for (unsigned int i = 0; i < s.size(); i++) {
        title << s[i] << " ";
      }
      title << endl;
    }
    
    oc.writeTitle(title.str());

    
        
    // loop over all trajectories 
    int frames=0;
    int time_counter=0;
    int init_counter=0;
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to   = args.upper_bound("traj");
    for (; iter != to; iter++) {
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()) {
	ic >> sys >> time;
	time_counter++;
	if (init_counter < first)
	  init_counter++;
	if (time_counter >= first && time_counter <= last ) {
	  oc.writeTimestep(time.steps(), time.time());
	  // only works for orthorombic boxes.check this
	  vec_box= (sys.box().K() + sys.box().L() + sys.box().M());
	  max =  atoms_gather.pos(index_center) + 0.5*vec_box;
	  min =  atoms_gather.pos(index_center) - 0.5*vec_box;
       
	  // loop over all selected atoms and gather

#ifdef OMP
#pragma omp parallel for  num_threads(threads)
#endif

	  for(int i = 0; i < atoms_gather.size(); i++) {
	    if (atoms_gather.pos(i)[0] > max[0]) 
	      atoms_gather.pos(i)[0]-=vec_box[0];
	    if (atoms_gather.pos(i)[1] > max[1])
	      atoms_gather.pos(i)[1]-=vec_box[1];
	    if (atoms_gather.pos(i)[2] > max[2])
	      atoms_gather.pos(i)[2]-=vec_box[2];
	    if (atoms_gather.pos(i)[0] < min[0])
	      atoms_gather.pos(i)[0]+=vec_box[0];
	    if (atoms_gather.pos(i)[1] < min[1])
	      atoms_gather.pos(i)[1]+=vec_box[1];
	    if (atoms_gather.pos(i)[2] < min[2])
	      atoms_gather.pos(i)[2]+=vec_box[2];
	  }
	  oc << atoms_gather;
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
	
//to do option to set pdb or cnf. 
// parallelize

