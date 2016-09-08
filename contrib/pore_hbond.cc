/**
 * @file pore_hbond.cc
 * Monitors the occurrence of hydrogen bonds inside a (0 0 0) centered pore aligned in the z axis
 */

#include <cassert>
#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/Hbond.h"
#include "../src/utils/groTime.h"
// added by jag
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/Box.h"
#include "../src/args/GatherParser.h"



using namespace gcore;
using namespace gio;
using namespace args;
using namespace utils;
using namespace std;
using namespace bound;
// added by jag
using namespace fit;

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "ref" << "DonorAtomsA" << "AcceptorAtomsA"
	 << "DonorAtomsB" << "AcceptorAtomsB" << "Hbparas"<< "time" 
	 << "massfile" << "traj"<<"first"<<"last"<<"pore"<<"doPore"
	 <<"min_z"<<"max_z"<<"atoms";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time          <time and dt>]\n";
  usage += "\t@DonorAtomsA    <atoms>\n";
  usage += "\t@AcceptorAtomsA <atoms>\n";
  usage += "\t@DonorAtomsB    <atoms>\n";
  usage += "\t@AcceptorAtomsB <atoms>\n";
  usage += "\t@Hbparas        <distance [nm],angle,min z, max z, pore radius**2; default: 0.25, 135,0,0,0>\n";
  usage += "\t[@massfile      <massfile>]\n";
  usage += "\t@traj           <trajectory files>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@pore    <pore atoms to be aligned with>\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  usage += "\t@doPore   <Define is hBonds are only compute for a 0 centered pore aligned in z axis>\n";
  usage += "\t@atoms  <atoms to follow>\n";


  try {
    Arguments args(argc, argv, knowns, usage);
    // Handle time
    Time time(args);
    // get the @first argument
    int first = args.getValue<int>("first");
    // get the @last argument
    int last = args.getValue<int>("last");
    // Read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    GromosForceField gff = it.forceField();
    // set the reference to fit to
    Reference reffit(&refSys);
    // read reference coordinates...
    if (args.count("ref") > 0) {
      InG96 ic(args["ref"]);
      ic.select("ALL");
      ic >> refSys;
      ic.close();
    } 
    else {
      InG96 ic;
      if (args.count("traj") > 0) {
	ic.open(args.lower_bound("traj")->second);
	ic.select("ALL");
	ic >> refSys;
	ic.close();
      } 
      else {
	throw gromos::Exception("pore_diffus", "no trajectory specified (@traj)");
      }
    }
    // The current system
    System sys(refSys);
    
    //fit atoms (the pore is aligned in z)
    AtomSpecifier fitatoms(refSys);
    {
      Arguments::const_iterator iter = args.lower_bound("pore");
      Arguments::const_iterator to = args.upper_bound("pore");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        fitatoms.addSpecifier(spec);
      }
    }
    reffit.addAtomSpecifier(fitatoms);
    RotationalFit rf(&reffit);


    // get pore atoms
    AtomSpecifier pore_atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("pore");
      Arguments::const_iterator to = args.upper_bound("pore");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        pore_atoms.addSpecifier(spec);
      }
    }
    if ( !pore_atoms.size() )
      throw gromos::Exception("pore_hbonds",
			      "No atoms to calculate pore dimensions!");
    
    //get atoms to count
    string atoms_count;
    AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
	atoms_count = spec;
        at.addSpecifier(spec);
      }
    }
    if ( !at.size() )
      throw gromos::Exception("pore_hbonds",
			      "No atoms to calculate!");

    

    
    //Get H-bonds parameters
    
    bool hbond3c = false;
    // get the paras
    vector<double> v_hbparas2c = args.getValues<double>("Hbparas", 5, false,
							Arguments::Default<double>() << 0.25 << 135.0 << 0.0 << 0.0 << 0.0);
    vector<double> v_hbparas3c;
    v_hbparas3c.resize(4);
    
    if (args.count("threecenter") >= 0) {
      hbond3c = true;
      v_hbparas3c = args.getValues<double>("threecenter", 4, false,
              Arguments::Default<double>() << 0.27 << 90.0 << 340.0 << 15.0);
    }
    HBPara2c hbparas2c = HB::mk_hb2c_paras_pore(v_hbparas2c);
    HBPara3c hbparas3c = HB::mk_hb3c_paras(v_hbparas3c);
    
    HB hb(sys, args, hbparas2c, hbparas3c);
    // initialize the calculation
    hb.init();

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
    
    InG96 ic;
    int frames=0;
    int time_counter=0;
    int init_counter=0;
    vector<int> loads_in_time;
    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj"); iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory

      while (!ic.eof()) {
        ic >> sys >> time;
	time_counter++;
	if (init_counter < first)
	  init_counter++;
	if (time_counter >= first && time_counter <= last ) {
	  frames++;
	  // do the fitting
	  (*pbc.*gathmethod)();
	  rf.fit(&sys);
	  static int frame = 1;
	  // get the number of atoms and break in case these numbers change from
	  // one frame to another
	  int numSoluAt = 0, numSolvAt = 0;
	  static int numSoluAt_old = -1, numSolvAt_old = -1;
	  int numSolu = sys.numMolecules();
	  int numSolv = sys.numSolvents();
	  for (int i = 0; i < numSolu; ++i) {
	    numSoluAt += sys.mol(i).numAtoms();
	  }
	  for (int i = 0; i < numSolv; ++i) {
	    numSolvAt += sys.sol(i).numAtoms();
	  }
	  if (numSoluAt_old != -1 && numSoluAt != numSoluAt_old) {
	    stringstream msg;
	    msg << "The number of solute atoms changed in " << iter->second.c_str() << ":\n"
		<< "             frame " << frame - 1 << ": " << numSoluAt_old << " solute atoms\n"
		<< "             frame " << frame << ": " << numSoluAt << " solute atoms\n"
		<< "       The calculation of hbond has been stopped therefore.";
	    throw gromos::Exception("hbond", msg.str());
	  }
	  if (numSolvAt_old != -1 && numSolvAt != numSolvAt_old) {
	    stringstream msg;
	    msg << "The number of solvent atoms changed in " << iter->second.c_str() << ":\n"
		<< "             frame " << frame - 1 << ": " << numSolvAt_old << " solvent atoms\n"
		<< "             frame " << frame << ": " << numSolvAt << " solvent atoms\n"
		<< "       The calculation of hbond has been stopped therefore.";
	    throw gromos::Exception("hbond", msg.str());
	  }
	  
	  numSoluAt_old = numSoluAt;
	  numSolvAt_old = numSolvAt;
	  
	  // get pore loading
	  int load_counter=0;
	  for(int i = 0; i < at.size(); i++) {
	    if (((at.pos(i)[2] >= (v_hbparas2c[2])) && (at.pos(i)[2] <= (v_hbparas2c[3]) ))  &&  
		( (at.pos(i)[0])*(at.pos(i)[0]) + (at.pos(i)[1])*(at.pos(i)[1]) <=  v_hbparas2c[4]) )
	      load_counter++;
	  }
	  loads_in_time.push_back(load_counter);

	  hb.settime(time.time());
	  hb.calc();
	  frame++;
	}	
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    // Final output of Data
    hb.printstatistics_pore(loads_in_time);
    
  } 
  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
