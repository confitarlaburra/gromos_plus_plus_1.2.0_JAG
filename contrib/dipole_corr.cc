#ifdef OMP
#include <omp.h>
#endif
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath> 
#include <algorithm> 

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Stat.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;
using namespace fit;

// JAG 11-02-14
//pore_P1: calculates P1 order parameter and  time autocorrelation function of particles located at the pore centre for 
// each individual dipole and a collective variabel defined as = Sum(dip_i) of all dipoles within the pore.
// It computes distributions along axis, full distributions, time-correlation functions and the time series of the collective
// dipole





//MAIN//

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "atoms"  
         << "traj"<<"first"
	 << "last"<< "pore"<<"window"<<"NatMol"
	 <<"threads";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@window    <time window>\n";
  usage += "\t@NatMol   <Number of atoms of molecule>\n";
  usage += "\t@threads    <thread number>\n";
  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);
       
    // get the @first argument
    int first = args.getValue<int>("first");
    
    // get the @last argument
    int last = args.getValue<int>("last");
    
    
     // get the @window argument
    int window = args.getValue<int>("window");

   
    // get NatmMol 
    int NatMol = args.getValue<int>("NatMol");

   
    
    // get the @threads argument
    int threads = args.getValue<int>("threads");


     //read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    {
      InG96 ic;
      if (args.count("traj") > 0) {
	ic.open(args.lower_bound("traj")->second);
	ic.select("ALL");
	ic >> refSys;
	ic.close();
      } else {
	throw gromos::Exception("dipole_corr", "no trajectory specified (@traj)");
      }
    }
    // The current system
    System sys(refSys);
    
                
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
      throw gromos::Exception("pore_diffus",
			      "No atoms to calculate adsorbed molecules!");
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
    
    vector<vector<Vec> > dipoles(at.size()); // 2D vector with all selected atoms and all postions of trajectory.
    vector<vector<Vec> > norm_dipoles(at.size());
    Vec z_axis (1.0,0.0,0.0);
    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    int frames=0;
    int time_counter=0;
    int init_counter=0;
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
	  frames++;
	  (*pbc.*gathmethod)();
	  //loop over pore molecules and get max and min z
	  for(int i = 0; i < (at.size()); i+=NatMol) {
	    Vec dipole(0.0,0.0,0.0);
	    for (int j =0; j < NatMol; j++) {
	      dipole+=at.charge(j+i)*at.pos(j+i);
	    }
	    dipole = dipole.normalize();
	    dipoles[i].push_back(dipole);
	    norm_dipoles[i].push_back(dipole.cross(z_axis).normalize());
	  }
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    


    //1d vector to auto-correlation 
    vector<double> P1_list(window,0.0);
    vector<double> norm_P1_list(window,0.0);
    // 1d vector to count for averaging
    vector<int> average_counter(window,0);
    
    
#ifdef OMP
#pragma omp parallel  num_threads(threads)
#endif
    for (int k =0; k < at.size(); k+=NatMol) { //loop trough molecules
      for(int it = 0; it < frames-window; it++) { // loop trough origins
	int window_counter = 0;
	for(int j = it; j < it+window; j++) { //loop through frames starting from it
	  P1_list[window_counter]+=(dipoles[k][j]).dot((dipoles[k][it]));
	  norm_P1_list[window_counter]+=(norm_dipoles[k][j]).dot((norm_dipoles[k][it]));
	  average_counter[window_counter]++;
	  window_counter++;
	}
      }
    }
    
    
    //output Dipole time-correlation data
    ofstream output;
    output.open("Mu_tcorr.dat");
    output<<"#<mu(t)dot mu(t0) of atoms "<<atoms_count<<endl;
    output<<"#Time[ps]          <mu*mu>          <norm_mu*norm_mu"<<endl;
    output.precision(4);
    for(int i = 0; i < P1_list.size(); i++) {
      output<<setw(15)<<i*time.dt()<<setw(15)<<(P1_list[i]/average_counter[i])
	    <<setw(15)<<(norm_P1_list[i]/average_counter[i])<<endl;
    }
    output.close();
    
    
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


// functions implementations//





