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


using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;


//JAG 24.02.2014
/* Program tha computes osmotic (pf) and diffusive permeabilities  of a single pore aligned in the z axis.
*/ 



// main

int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "pore" << "atoms"  
         << "traj" << "ref" << "offset" << "pore_radius"
	 << "threads"<< "first" <<"last"<<"window";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@threads   <number of threads  >\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  usage += "\t@window   <time window for average [ps]  >\n";
  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    // get the @offset argument
        
    // get the @offset argument
    double offset=0.0;
    if ( args.count("offset") > 0)
      offset=args.getValue<double>("offset");
    // get the @first argument
    int first = args.getValue<int>("first");
    
    // get the @last argument
    int last = args.getValue<int>("last");

    // get the @pore_radius argument
    double pore_radius = args.getValue<double>("pore_radius");
    //get threads
    int threads;
    if (args.count("threads") > 0) 
      threads=args.getValue<int>("threads");
    else
      threads=1;      

    int window;
    if (args.count("window") > 0) 
      window=args.getValue<int>("window");
    else
      window=100;      


    //read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    // set the reference to fit to
    Reference reffit(&refSys);
    
    // read reference coordinates...
    if (args.count("ref") > 0) {
      InG96 ic(args["ref"]);
      ic.select("ALL");
      ic >> refSys;
      ic.close();
      
    } else {
      InG96 ic;
      if (args.count("traj") > 0) {
	ic.open(args.lower_bound("traj")->second);
	ic.select("ALL");
	ic >> refSys;
	ic.close();
      } else {
	throw gromos::Exception("frameout", "no trajectory specified (@traj)");
      }
    }

   
    // The current system
    System sys(refSys);
    
    
    //fit atoms
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
      throw gromos::Exception("pore_ads",
			      "No atoms to calculate pore dimensions!");
    
        
    //get atoms to count
    string atoms_count;
    AtomSpecifier count_atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
	atoms_count = spec;
        count_atoms.addSpecifier(spec);
      }
    }
    if ( !count_atoms.size() )
      throw gromos::Exception("pore_ads",
			      "No atoms to calculate adsorbed molecules!");
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
    
    // 1d vector to store actual and old z positions of atoms to count in the pore
    vector<double> z_list(count_atoms.size(),0), z_list_old(count_atoms.size(),0);
    //vector<double> z_list_old(count_atoms.size(),0);
    vector<double> x_list(count_atoms.size(),0),x_list_old(count_atoms.size(),0);
    //vector<double> x_list_old(count_atoms.size(),0);
    vector<double> y_list(count_atoms.size(),0),y_list_old(count_atoms.size(),0);
    // list to store dn values
    vector<double> dn_list;
    //define rest of the variables
   
    double max_z,min_z,center_x,center_y,rad2,rad;
    double dzAll=0.0;
    double dz=0.0;
    double dn=0.0;
    double n =0.0;
    double n_old=0.0;
    int pLoading=0;
    int up=0;
    int down=0;
    ofstream output;
    output.open("Nt.out");
    output << "# Collective variable n  and permeation events of "<<atoms_count<<endl;
    output << "# Time[ps]           [n]                Load               z+              z-            Total"<<endl;
    output.precision(3);
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
	  // do the fitting
	  (*pbc.*gathmethod)();
	  rf.fit(&sys);
	  //loop over pore molecules and get max and min z
	
	  min_z = pore_atoms.pos(0)[2];
	  max_z = pore_atoms.pos(0)[2];
	  		
	  	
	  // Loop trough all pore atoms to get min and max z of the pore
#ifdef OMP
#pragma omp parallel for  num_threads(threads)	\
  reduction(+: center_x,center_y)
#endif
	  for(int i = 0; i < pore_atoms.size(); i++) {
	    center_x+=pore_atoms.pos(i)[0];
	    center_y+=pore_atoms.pos(i)[1];
	    if(pore_atoms.pos(i)[2] > max_z)
	      max_z = pore_atoms.pos(i)[2];
	    if(pore_atoms.pos(i)[2] < min_z)
	      min_z = pore_atoms.pos(i)[2];
	  }
	  // set COG in the x-y plane
	  center_x /= pore_atoms.size();
	  center_y /= pore_atoms.size();
	  // set maximal radial distance in x-y plane
	  rad2= pore_radius*pore_radius;
	  	
	  // loop through selected atoms and collect data

	  
	  
	  z_list_old=z_list;
	  x_list_old=x_list;
	  y_list_old=y_list;
	  //cout<<"frames "<<frames<<endl;
	  
#ifdef OMP
#pragma omp parallel for  num_threads(threads)		
#endif
	  for(int i = 0; i < count_atoms.size(); i++) {
	    z_list[i]=count_atoms.pos(i)[2];
	    x_list[i]=count_atoms.pos(i)[0]- center_x;
	    y_list[i]=count_atoms.pos(i)[1]- center_y;
	  }
	  //counters set to zero
	  dzAll=0.0;
	  pLoading =0;
	  dn=0.0;
	  /*iterates over each water molecule, locates the ones that are and were inside the pore
	    and computes their unidimensional (z) displacement, dz
	    dz is accumulated in dzAll
	    also counts the amount of waters molecules inside the pore*/

#ifdef OMP
#pragma omp parallel for  num_threads(threads)	\
  reduction(+: dzAll,pLoading)
#endif
	  
	  for(int i = 0; i < count_atoms.size(); i++) 
	    if ( (z_list[i] >= (min_z - offset)) && z_list[i] <= (max_z + offset) )  
	      if ((z_list_old[i] >= (min_z - offset)) && z_list_old[i] <= (max_z + offset) ) { 
		double rad2_actual= x_list[i]*x_list[i] + y_list[i]*y_list[i];
		double rad2_old= x_list_old[i]*x_list_old[i] + y_list_old[i]*y_list_old[i];
		// select molecules within a max radius
		if( rad2_actual <=  rad2) 
		  if( rad2_old  <=  rad2)  {
		    dz =z_list[i]-z_list_old[i];
		    dzAll+=dz;
		    pLoading++;
		  }
	      }
	  
	  
	  //computes collective variable dn and integrates it (summing up) getting n(t)
	  if (frames==1)
	    dn=0.0;
	  else
	    dn = dzAll/((max_z + offset)- (min_z - offset));
	  if (frames==1)
	    dn_list.push_back(0.0);
	  else
	    dn_list.push_back(dn);
	  n += dn;
	  //Get total permeation events in the -z direction
	  if ((n-n_old) <= -1) {
	    n_old=n;
	    down++;
	  }
	  //Get total permeation events in the +z direction
	  if ((n-n_old) >= 1) {
	    n_old=n;
	    up++;
	  }
	  output << time;
	  output<<setw(15)<<n<<setw(15)<<pLoading<<setw(15)<<up<<setw(15)<<down<<setw(15)<<up+down<<endl;
	  
	}
      }
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    
    output.close();
    
    //1d vector to store <N2>
    vector<double> N2_list(window,0.0);
    // 1d vector to count for averaging
    vector<int> average_counter(window,0);
    
    //double loop for windowed multiple origins

#ifdef OMP
#pragma omp parallel for  num_threads(threads)	
#endif
    for(int it = 0; it < frames-window; it++) {
      int window_counter = 0;
      double n_loop=0;
      double N2=0;
      for(int j = it; j < it+window; j++) {
	n_loop+=dn_list[j];
	//computes the msd of n, N2 and accumualtes for every timewindow frame
	N2=n_loop*n_loop;
	N2_list[window_counter]+=N2;
	average_counter[window_counter]++;
        window_counter++;
      }
    }
   

    output.open("N2_msd.dat");
    output<<"#Time[ps]   <N2>"<<std::endl;
    output.precision(3);
    double t=0.0;
    output<<setw(15)<<t<<setw(15)<<t<<endl;
    for(int i = 0; i < N2_list.size(); i++) {
      output<<setw(15)<<i+1<<setw(15)<<(N2_list[i]/average_counter[i])<<endl;
    }
    output.close();

        
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

//functions and objects implementations

