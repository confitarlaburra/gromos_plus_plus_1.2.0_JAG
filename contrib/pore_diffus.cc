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
//pore_diffus: calculates diffusion and mean residence time inside a cylindrical pore.
// it only considers the closest particle to the center of the pore.
// It needs a reference structure with the pore aligned in z. 


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
  knowns << "topo" << "pbc" << "time" << "dim" << "atoms"  
         << "traj"<<"ref"<<"offset"<<"pore_radius"<<"first"
	 << "last"<< "pore"<<"window"<<"origin"<<"WriteTraj"
	 <<"filtered"<<"threads";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@dim    <dimensions to consider>\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@window    <time window>\n";
  usage += "\t@filtered    <only consider a molecule at center once if it did not lefet the pore between measurements for residence times>\n";
  usage += "\t@WriteTraj    <write trajectories of particles starting from pore center: 0:no, 1:yes>\n";
  usage += "\t@threads    <thread number>\n";
  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    // get the relevant dimensions
    int ndim = 3;
    int dim[3] = {0, 1, 2};
    {
      Arguments::const_iterator iter = args.lower_bound("dim");
      if (iter != args.upper_bound("dim")) {
        ndim = 0;

        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
        iter++;
      }
      if (iter != args.upper_bound("dim")) {
        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
        iter++;
      }
      if (iter != args.upper_bound("dim")) {
        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
      }
    }


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
    int window = args.getValue<int>("window");

    // get the @window argument
    int WriteTraj = args.getValue<bool>("WriteTraj");
   
     // get the @filtered argument
    int filtered = args.getValue<bool>("filtered");
   
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
      throw gromos::Exception("pore_diffus",
			      "No atoms to calculate adsorbed molecules!");
    
      AtomSpecifier old_at = at;
      old_at.setSystem(oldsys);

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
   


    
    vector<vector<Vec> > position(at.size()); // 2D vector with all selected atoms and all postions of trajectory.
    vector<double> maxs_z; // vector with  min and max of pore axial axis
    vector<double> mins_z;
    vector<double> centers_x;
    vector<double> centers_y;
    vector<int> atoms_at_center; // vector with indexes of closest atoms to the pore center
    vector<vector <bool> > in_out(at.size()); // 2D vector with all selected atoms tracking if atoms were inside (1) or outside the pore (0)
    vector<vector<int > >loads_in_time;
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
	  
	  mins_z.push_back(min_z);
	  maxs_z.push_back(max_z);
	  

	  // set COG in the x-y plane
	  center_x /= pore_atoms.size();
	  center_y /= pore_atoms.size();
	
	  centers_x.push_back(center_x);
	  centers_y.push_back(center_y);
	  
	  // loop through selected atoms
	  vector<int> load_indexes;
	  for(int i = 0; i < at.size(); i++) {
	    at.pos(i) = pbc->nearestImage(old_at.pos(i), at.pos(i), sys.box());
	    position[ i ].push_back(at.pos(i));
	    old_at.pos(i) = at.pos(i);
	    if (((at.pos(i)[2] >= (min_z - offset)) && (at.pos(i)[2] <= (max_z + offset) ))  &&  
	      		( (at.pos(i)[0]-center_x)*(at.pos(i)[0]-center_x) + (at.pos(i)[1]-center_y)*(at.pos(i)[1]-center_y) <=  rad2) ) {
	      load_indexes.push_back(i);
	      in_out[i].push_back(true);
	    } else
	      in_out[i].push_back(false);
	  }
	  
	  loads_in_time.push_back(load_indexes);
	  if (load_indexes.size() > 0) 
	    atoms_at_center.push_back(load_indexes[find_center(load_indexes,at,min_z,max_z)]);
	  else
	    atoms_at_center.push_back(-1);
	  //cout << time;
	  //cout<<setw(9)<<load_indexes.size()<<" "<<load_indexes[find_center(load_indexes,at,min_z,max_z)]<<endl;
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }


	  
    
    //1d vector to store msd
    vector<double> MSD_list(window,0.0);
    // 1d vector to count for averaging
    vector<int> average_counter(window,0);
    // 1d vector to count residence time for averaging
    vector<int> residence_times;
    cout<<"#\tNow computing ...."<<endl;
    ofstream output;
  
    vector<int> filtered_indexes = filter_array(atoms_at_center, in_out); 
    //for(int i = 0; i < filtered_indexes.size(); i++) {
    //  cout<<" "<<atoms_at_center[i]<<" "<<filtered_indexes[i]<<endl;
    //}
 

    
    cout<<"# Frames: "<<frames<<endl;
    // mean residence time code
    //int filtered =1;
    for(int i = 0; i < frames;  i++) {
      int index =0;
      if (filtered)
	index = filtered_indexes[i];
      else
	index = atoms_at_center[i];
      int residence_time=0.0;
      int test;
      if (index != -1) {
	
	if (WriteTraj) {
	  string name1 = at.toString(index);
	  string name2 = static_cast<ostringstream*>( &(ostringstream() << i) )->str(); //convert int into string
	  string name3 = name2 + "." + name1 + ".out";
	  output.open(name3.c_str());
	  output.precision(3);
	  output<<"#Trajectory from center of atom: "<<at.toString(index)<<endl;
	  output<<"#Time           z[nm]"<<endl;
	  output<<setw(10)<<i<<setw(10)<<(position[index][i])[2]<<endl;
	}
	for(int it = i+1; it < frames; it++) {
	  if ( ( (position[index][it])[2] >= (mins_z[it] - offset) ) && ( (position[index][it])[2] <= (maxs_z[it] + offset) ) ) {
	    if (WriteTraj) 
	      output<<setw(10)<<it<<setw(10)<<(position[index][it])[2]<<endl;
	  }
	  else {
	    residence_times.push_back(it-i);
	    cout<<"\natom "<<at.toString(index)<<" at center in frame "<<i<<" left at "<<it<<endl;
	    if (WriteTraj) 
	      output<<setw(10)<<it<<setw(10)<<(position[index][it])[2]<<endl;
	    break;
	  }
	}
	if (WriteTraj)
	  output.close();
      }
    }
    

   
    
    output.open("residence_times.dat");
    output<<"#Residence times [ps] of atoms "<<atoms_count<<endl;
    output.precision(3);
    
    Stat<double> RT_error;
    if ( !residence_times.size() ) {
      cout<<"\tNo atoms of "<<atoms_count<<" left the pore from its center in "<<frames<<" snapshots"<<endl;
      output<<"\tNo atoms of "<<atoms_count<<" left the pore from its center in "<<frames<<" snapshots"<<endl;
    }
    else  {
      for(int i = 0; i < residence_times.size(); i++) {
	output<<setw(15)<<residence_times[i]<<endl;
	RT_error.addval(residence_times[i]);
      }
      output<<"#Average = "<<RT_error.ave()<<endl;
      output<<"#RMSD    = "<<RT_error.rmsd()<<endl;
      output<<"#Error   = "<<RT_error.ee()<<endl;
    }
    
    output.close();

   
     //Windowed MSD cycles with multiple origins 

#ifdef OMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for(int it = 0; it < frames-window; it++) {
      int window_counter = 0;
      double square_d = 0;
      for(int j = it; j < it+window; j++) { //loop through frames starting from it
	for (int h =0; h < loads_in_time[it].size() ; h++) { // loop through loads at origin time it
	  int i=loads_in_time[it][h];  // select particle i loaded at time origin it
	  square_d = 0;
	  if (in_out[i][j]) { // checdk if it was loaded at time j
	    for(int k = 0; k < ndim; k++) {
	      const double d = (position[i][j])[dim[k]]-(position[i][it])[dim[k]];
	      square_d += d * d;
	    }
	    MSD_list[window_counter]+=square_d;
	    average_counter[window_counter]++;
	  }
	}
	window_counter++;
      }
    }


    
    //output MSD data
    
    output.open("msd.dat");
    output<<"#MSD of atoms "<<atoms_count<<endl;
    output<<"#Time[ps]          <MSD>"<<endl;
    output.precision(3);
    double t=0.0;
    output<<setw(15)<<t<<setw(15)<<t<<endl;
    
    for(int i = 0; i < MSD_list.size(); i++) {
      output<<setw(15)<<i*time.dt()<<setw(15)<<(MSD_list[i]/average_counter[i])<<endl;
    }
    output.close();      




    //only use COG of tube a orgin for msd calculation (less data)
 
    /*
    for(int it = 0; it < frames-window; it++) {
      int window_counter = 0;
      double square_d = 0;
      int i = atoms_at_center[it]; //index of atom at centre of pore at that particular snapshot
      if (i!=-1)
	for(int j = it; j < it+window; j++) { //loop through frames starting from it
	  square_d = 0;
	  if ( ( (position[i][j])[2] >= (mins_z[j] - offset) ) && ( (position[i][j])[2] <= (maxs_z[j] + offset) )  
	       && ( (position[i][it])[2] >= (mins_z[it] - offset) ) && ( (position[i][it])[2] <= (maxs_z[it] + offset) )   ) {
	    for(int k = 0; k < ndim; k++) {
	      const double d = (position[i][j])[dim[k]]-(position[i][it])[dim[k]];
	      square_d += d * d;
	    }
	    MSD_list[window_counter]+=square_d;
	    average_counter[window_counter]++;
	    window_counter++;
	  }
	  else
	    break;
	}
    }
    
    //full cyle considering all atoms
	
    for(int it = 0; it < frames-window; it++) {
      int window_counter = 0;
      double square_d = 0;
      //cout<<"frame "<<it<<endl;
      for(int j = it; j < it+window; j++) { //loop through frames starting from it
	//cout<<" subframe "<<j<<" window "<<window_counter<<endl;
	int loads =0;
	for (int i =0; i < at.size(); i++) {
	  square_d = 0;
	  if ( ( (position[i][it])[2] >= (mins_z[it] - offset) ) && ( (position[i][it])[2] <= (maxs_z[it] + offset) ) ) {
	    if (  ( (position[i][it])[0]-centers_x[it])*( (position[i][it])[0]-centers_x[it])  +  ( (position[i][it])[1]-centers_y[it])*( (position[i][it])[1]-centers_y[it]) <= rad2) {
	      if ( ( (position[i][j])[2] >= (mins_z[j] - offset) ) && ( (position[i][j])[2] <= (maxs_z[j] + offset) ) ) {
		if (  ( (position[i][j])[0]-centers_x[j])*( (position[i][j])[0]-centers_x[j])  +  ( (position[i][j])[1]-centers_y[j])*( (position[i][j])[1]-centers_y[j]) <= rad2) {
		  //loads++;
		  //cout << mins_z[j]<<" "<< maxs_z[j]<<" "<<endl;
		  for(int k = 0; k < ndim; k++) {
		    const double d = (position[i][j])[dim[k]]-(position[i][it])[dim[k]];
		    square_d += d * d;
		  }
		  MSD_list[window_counter]+=square_d;
		  average_counter[window_counter]++;
		}
	      }
	    }
	  }
	}
	//cout<<"loads "<<loads<<endl;
	window_counter++;
      }
    }
    */

    
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


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
	    tube_state=in_out[element][k]; //go back in time and check if in between j and i it left the pore (tube_state==0) 
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

