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

// JAG 28-02-14
//pore_dist_corr: Computes displacement correlation of all atoms pairs loaded in pore
// For more details please see J.Chem.Phys 126, 124704 (2007)


//MAIN//

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time"<< "atoms"  
         << "traj"<<"ref"<<"offset"<<"pore_radius"<<"first"
	 << "last"<< "pore"<<"origin"<<"bins"<<"full"<<"threads";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@bins    <number of bins for distances correlation\n";
  usage += "\t@full    <compute  the distance correlation for all space";
  usage += "\t@threads    <thread number >";
  
  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

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
   
     // get the @window argument
    //int window = args.getValue<int>("window");

    // get the @window argument
    // int filtered = args.getValue<int>("filtered");

    // get the @bins argument
    int bins = args.getValue<int>("bins");
        
    // get the @bins argument
    bool full = args.getValue<bool>("full");

    // get the @threads argument
    int threads = args.getValue<int>("threads");

    
    

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
	throw gromos::Exception("pore_diffus", "no trajectory specified (@traj)");
      }
    }

   
    // we always need the old coordinates to take the nearest image
    System oldsys(refSys);
    
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
      throw gromos::Exception("pore_dist_corr",
			      "No atoms to calculate!");

    // for ease of looping, we make three of these atomspecifiers
    // one for each system
    AtomSpecifier ref_at = at;
    ref_at.setSystem(refSys);
    AtomSpecifier old_at = at;
    old_at.setSystem(oldsys);
    
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
   


    
    vector<vector<Vec> > position(at.size()); // 2D vector with all selected atoms and all postions of trajectory.
    //vector<vector<Vec> > displacements(at.size()); // 2D vector with all selected atoms and all displacements of trajectory.
    vector<vector<int > >loads_in_time;
    // 2D vector with all selected atoms tracking if atoms were inside (true) or outside the pore (false)
    vector<vector <bool> > in_out(at.size()); 
    double  max_z,min_z,center_x,center_y;
    double rad2= pore_radius*pore_radius;
    int pore_load =0;

        
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
	  
	  // loop through selected atoms
	  vector<int> load_indexes;
	  for(int i = 0; i < at.size(); i++) {
	    at.pos(i) = pbc->nearestImage(old_at.pos(i), at.pos(i), sys.box());
	    position[ i ].push_back(at.pos(i));
	    if (((at.pos(i)[2] >= (min_z - offset)) && (at.pos(i)[2] <= (max_z + offset) ))  &&  
		( (at.pos(i)[0]-center_x)*(at.pos(i)[0]-center_x) + (at.pos(i)[1]-center_y)*(at.pos(i)[1]-center_y) <=  rad2) ) {
	      load_indexes.push_back(i);
	      in_out[i].push_back(true);
	    } else
	      in_out[i].push_back(false); 
	  }
	  loads_in_time.push_back(load_indexes);
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }


	  
    
    //1d vector to store distance correlation
    
    vector<float> corr_dist(bins);;
    // 1d vector to count for averaging
    vector<int> average_counter(bins);
    // init variable
    float bin_size =(abs(max_z-min_z))/bins;
    cout<<"#\tNow computing ...."<<endl;
    cout<<"# Frames: "<<frames<<endl;
    // multiple origins
    if (!full)

#ifdef OMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
      for(int it = 0; it < frames; it++) {
	for(int j = it; j < frames; j++) { //loop through frames starting from it. i.e multiple origins
	  // Double loop for displacements of pairs i and m pairs within the the tube
	  for (int h =0; h < loads_in_time[it].size() ; h++) {
	    int i =loads_in_time[it][h]; // index of atom h loaded in tube at time t
	    //Check that particle i  is within tube at frame  j
	    if (in_out[i][j]) {
	      float disp_i =(position[i][j])[2] - (position[i][j-1])[2]; // displacement for particle i in t+dt
	      for (int k =h+1; k < loads_in_time[it].size() ; k++) {
		int m =loads_in_time[it][k]; // index of atom k loaded in tube at time it
		//Check that particle m  is within tube at frame j
		if (in_out[m][j]) {
		  float disp_m=(position[m][j])[2]-(position[m][j-1])[2];
		  float axial_distance =abs( (position[m][j])[2] - (position[i][j])[2]);
		  float step= floor(axial_distance/bin_size);
		  int bin=int(step);
		  if (bin >= 0 && bin < bins) { 
		    corr_dist[bin]+=disp_i*disp_m; 
		    average_counter[bin]++;
		  }
		}
	      }
	    }
	  }
	}
      }


    

    // full correlation code 
    if (full)
      bin_size = 4.8/bins;
    if (full)      
#ifdef OMP
#pragma omp parallel for  schedule(dynamic) num_threads(threads)
#endif
      for(int it = 0; it < frames; it++) {
	cout<<"outer "<<it<<endl;
	for(int j = it; j < frames; j++) {
	  for (int i =0; i < at.size() ; i++) {
	    float disp_i =(position[i][j])[2] - (position[i][j-1])[2];
	    for (int k=i+1; k<at.size(); k++) {
	      float disp_k=(position[k][j])[2]-(position[k][j-1])[2];
	      float axial_distance =abs( (position[k][j])[2] - (position[i][j])[2]);
	      float step= floor(axial_distance/bin_size);
	      int bin=int(step);
	      if (bin >= 0 && bin < bins) { 
		corr_dist[bin]+=disp_i*disp_k; 
		average_counter[bin]++;
	      }
	    }
	  }
	}
      }
       
    //output correlation data
    ofstream output;   
    output.open("disp_corr.dat");
    output<<"#Displacement correlation  of atoms "<<atoms_count<<endl;
    output<<"#Axial Distance          <drj*dri>"<<endl;
    output.precision(3);
    for(int i = 0; i < corr_dist.size(); i++) {
      output<<setw(15)<<i*bin_size + bin_size*0.5<<setw(15)<<(corr_dist[i]/average_counter[i])<<endl;
    }
    output.close();      
    
     
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



