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

// JAG 11-04-14
//pore_surv_P: calculates survival probabilities employing Q(t)=sum<pit(Pi)> where
// Pi is the dirac delta function, 1 inside pore 0 elsewhere.
// For more details please see J.Chem.Phys 126, 124704 (2007)


// Functions declarations
/*
  Finds the closest molecule to the center of the pore in the z axis.
  input:
  indexes array = array with molecules at pore center at each frame.
  atoms = selected atoms
  min_z, max_x = min and max coordinates   
*/

int find_center( vector<int> &indexes, AtomSpecifier &atoms, double &min_z, double &max_z);

/*
  Filter arrays for mean residence time calculations from pore center.considers a molecule at the center only once if
  it did not left pore between the measurements
  input:
  indexes array = array with molecules at pore center at each frame.
  in_out= 2D array which contails all selected molecules and for each a counter indicating if they are at the pore (1)
  or outside the pore (0) for each frame of the simulation.
*/

vector<int> filter_array(vector<int> &indexes , vector < vector<bool> > &in_out) ;





//MAIN//

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time"<< "atoms"  
         << "traj"<<"ref"<<"offset"<<"pore_radius"
	 <<"first"<< "last"<< "pore"<<"threads"
	 <<"full"<<"filtered";

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
  usage += "\t@threads    <thread number >";
  usage += "\t@full    <compute survival times of all loads or just molecules at center>";
  usage += "\t@filtered    <use filtered loads for molecules at center (consider only once unless it lef the pore)>";
  
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
   
    // get the @threads argument
    int threads = args.getValue<int>("threads");

    // get the @full argument
    bool full = args.getValue<bool>("full");
    
    // get the @filtered argument
    bool filtered = args.getValue<bool>("filtered");
    
    
    

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
    vector<vector<double> > Avg_Pi(at.size()); //2D vector with all selected atoms and average delta functions per time.
    vector<vector<double> > Avg_counter(at.size()); // Counter to get averages
    vector<vector<int > >loads_in_time;
    vector<int> atoms_at_center; // vector with indexes of closest atoms to the pore center
    // 2D vector with all selected atoms tracking if atoms were inside (true) or outside the pore (false)
    vector<vector <bool> > in_out(at.size()); 
    double  max_z,min_z,center_x,center_y;
    double rad2= pore_radius*pore_radius;
    
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
	    // Just to fill vectors
	    Avg_Pi[i].push_back(0); 
	    Avg_counter[i].push_back(0);
	    // if inside pore
	    if (((at.pos(i)[2] >= (min_z - offset)) && (at.pos(i)[2] <= (max_z + offset) ))  &&  
		( (at.pos(i)[0]-center_x)*(at.pos(i)[0]-center_x) + (at.pos(i)[1]-center_y)*(at.pos(i)[1]-center_y) <=  rad2) ) {
	      load_indexes.push_back(i);
	      in_out[i].push_back(true);
	    } else
	      in_out[i].push_back(false);
	  }
	  loads_in_time.push_back(load_indexes);
	  if (load_indexes.size() > 0 && !full)
	    // Get particle closest to the axial center of the pore
	    atoms_at_center.push_back(load_indexes[find_center(load_indexes,at,min_z,max_z)]);
	  else
	    atoms_at_center.push_back(-1);
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }

        
    cout<<"#\tNow computing ...."<<endl;
    cout<<"# Frames: "<<frames<<endl;
    
    // consider all atoms loaded in the tube
    if (full)
#ifdef OMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
      for(int it = 0; it < frames; it++) {
	int time=0;
	for(int j = it; j < frames; j++) { //loop through frames starting from it. i.e multiple origins
	  for (int h =0; h < loads_in_time[it].size() ; h++) {
	    int i =loads_in_time[it][h]; // index of atom h loaded in tube at time t
	    if (in_out[i][j])  {//Check that particle i  is within tube at frame j
	      Avg_Pi[i][time]+=1;
	      // Avg_counter[i][time]++;
	    }
	  }
	  time++;
	}
      }

    // consider only atoms at the centre of the pore
    if (!full) {
      vector<int> filtered_indexes;
      // only consider particles located at pore center once unless it has left the pore 
      if (filtered)
	filtered_indexes = filter_array(atoms_at_center, in_out);
      else
	filtered_indexes=atoms_at_center;
#ifdef OMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
      for(int it = 0; it < frames; it++) {
	int time=0;
	int i =filtered_indexes[it]; // index of atom h at center of  tube at time it
	if (i!=-1)
	  for(int j = it; j < frames; j++) { //loop through frames starting from it. i.e multiple origins
	    if (in_out[i][j]) { //Check that particle i  is within tube starting at frame j 
	      Avg_Pi[i][time]+=1;
	      //Avg_counter[i][time]++;
	    }
	    time++;
	  }
      }
    }
    

    // Get Survival probability Q(t)   
    vector<double> Qt(frames);
    for(int j = 0; j < at.size(); j++) {
      for (int i=0; i<frames; i++) {
	//if (Avg_counter[j][i]>0 && Avg_counter[j][i] > 0 )
	  //Qt[i]+=(Avg_Pi[j][i])/(Avg_counter[j][i]);
	  Qt[i]+=(Avg_Pi[j][i]);
      }
    }
   
    // write output 
    ofstream output;   
    output.open("Surv_prob.dat");
    output<<"#Survival probability of atoms "<<atoms_count<<endl;
    output<<"#Time                Q(t)"<<endl;
    for (int i=0; i<frames; i++) 
      output<<setw(10)<<i<<setw(15)<<Qt[i]/Qt[0]<<endl;
    output.close();
    
    
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

//END//
// functions implementations//


// Finds the closest molecule to the center of the pore in the z axis.

int find_center( vector<int> &indexes, AtomSpecifier &atoms, double &min_z, double &max_z) {
  double z  = 0.0;
  double center_z = (max_z + min_z)*0.5;
  double small = abs ( atoms.pos(indexes[0])[2]-center_z);
  int index = 0;
  for (int i = 0; i < indexes.size(); i++) {
    z=abs(atoms.pos(indexes[i])[2]-center_z);
    if (z<small) {
      small=z;
      index =i;
    }
  }
  return index;
}

//Filter arrays for mean residence time calculations from pore center
//input:
// indexes array = array with molecules at pore center at each frame.
// in_out= 2D array which contains all selected molecules and for each, a counter indicating if they are at the pore (1)
// or outside the pore (0) for each frame of the simulation.
 
vector<int> filter_array(vector<int> &indexes , vector < vector<bool> > &in_out) {
  vector<int> filtered_array(indexes.size(),-1);
  vector<int> repeats(in_out.size());
  int element;
  int tube_state;
  int counter =0;
  for (int i=0; i <indexes.size(); i++) {
    if (indexes[i] != -1) { 
      element =indexes[i];
      repeats[element]++;
      if (repeats[element] < 2)
	filtered_array[i]=element; //always add element when observed for the first time 
      for (int j=i+1; j <indexes.size(); j++) {
	if (element == indexes[j]) {  //if it appears again in the list
	  for (int k=j; k >=i; k--) { 
	    tube_state=in_out[element][k]; //go back in time and check if in between j and i it left the pore (tube_state==false) 
	    if (tube_state == false)
	      counter++;      
	  }
	  if (counter > 0) 
	    filtered_array[j]=element; // if the former occurred at least once, also add these repeated element 
	}
      }
    }
  }
  return filtered_array;
}




