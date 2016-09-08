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
         << "traj"<<"ref"<<"offset"<<"pore_radius"<<"first"
	 << "last"<< "pore"<<"window"<< "bins"<<"NatMol"<<"WriteDistBin"
	 << "FullCorrelation"<<"threads";

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
  usage += "\t@bins    <number of bins for axial distribution>\n";
  usage += "\t@window    <time window>\n";
  usage += "\t@NatMol   <Number of atoms of molecule>\n";
  usage += "\t@WriteDistBin   <Write P1 distributions for each bin of the pore>\n";
  usage += "\t@FullCorrelation   <Compute mu autocorrelation function for all molecules in the box>\n";
  usage += "\t@threads    <thread number>\n";
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
    int window = args.getValue<int>("window");

    // get bins
    int bins = args.getValue<int>("bins");

    // get NatmMol 
    int NatMol = args.getValue<int>("NatMol");

    // get WriteDistBin 
    int WriteDistBin = args.getValue<int>("WriteDistBin");

    //get FullCorrelation
    int  FullCorrelation = args.getValue<int>("FullCorrelation");

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
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
   


    
    vector<vector<Vec> > dipoles(at.size()); // 2D vector with all selected atoms and all postions of trajectory.
    vector<vector<Vec> > norm_dipoles(at.size());// vector normal to dipole
    vector<Vec> Total_dipoles;
    vector<Vec> norm_Total_dipoles; // vector normal to total dipole
    //vector<vector<Vec> > position(at.size()); // 2D vector with all selected atoms and all dipoles  of trajectory.
    vector<double> maxs_z; // vector with  min and max of pore axial axis
    vector<double> mins_z;
    vector<double> centers_x;
    vector<double> centers_y;
    vector<vector<int > >loads_in_time;
    vector<vector <bool> > in_out(at.size()); // 2D vector with all selected atoms tracking if atoms were inside (1) or outside the pore (0)
    // 1d vector to store histogram  along pore axis
    vector<vector<double > > hists1D(bins); 
    vector<vector<double > > norm_hists1D(bins);
    double  max_z,min_z,center_x,center_y;
    double rad2= pore_radius*pore_radius;
    double bin_size;
    Distribution P1(-1.0,1.0,bins);
    Distribution P1_collective(-1.0,1.0,bins);
    
    Distribution norm_P1(-1.0,1.0,bins);
    Distribution norm_P1_collective(-1.0,1.0,bins);
    
    Vec z_axis (0.0,0.0,1.0);
    Vec x_axis (1.0,0.0,0.0);
    int step,bin;
    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    int frames=0;
    int time_counter=0;
    int init_counter=0;
    ofstream output("P1_collective_TS.dat");
    output<<"#Collective dipole u norm trajectory and collective P1 of "<<atoms_count<<endl;
    output<<"#Time                  u abs        P1 collec     P1(x) collec  norm"<<endl;    
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
	  //if (args.count("ref") > 0)
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
	  
	  bin_size =((max_z+offset)-(min_z-offset))/(bins);

	  // set COG in the x-y plane
	  center_x /= pore_atoms.size();
	  center_y /= pore_atoms.size();
			
	  centers_x.push_back(center_x);
	  centers_y.push_back(center_y);

	  // loop through selected atoms
	  vector<int> load_indexes;
	  Vec Total_dipole(0.0,0.0,0.0);
	  Vec norm_Total_dipole(0.0,0.0,0.0);
	  double dot_collective=0.0;
	  for(int i = 0; i < (at.size()); i+=NatMol) {
	    Vec dipole(0.0,0.0,0.0);
	    Vec dipole2(0.0,0.0,0.0);
	    for (int j =0; j < NatMol; j++) {
	      dipole+=at.charge(j+i)*at.pos(j+i);
	      dipole2+=at.charge(j+i)*at.pos(j+i);
	    }
	    //Total_dipole+=dipole2; // if the total dipole of the box needs to be computed
	    dipole = dipole.normalize();
	    dipoles[i].push_back(dipole);
	    norm_dipoles[i].push_back(dipole.cross(z_axis).normalize());
	    //position[i].push_back(at.pos(i));
	    if (((at.pos(i)[2] >= (min_z - offset)) && (at.pos(i)[2] <= (max_z + offset) )) &&  
		( (at.pos(i)[0]-center_x)*(at.pos(i)[0]-center_x) + (at.pos(i)[1]-center_y)*(at.pos(i)[1]-center_y) ) <=  rad2 ) {
	      load_indexes.push_back(i);
	      Total_dipole+=dipole2;
	      in_out[i].push_back(true);
	      step= floor( (at.pos(i)[2]-(min_z-offset))/(bin_size));
	      bin=int(step);
	      if (bin >= 0 && bin < bins) {
		hists1D[bin].push_back(dipole.dot(z_axis));
		norm_hists1D[bin].push_back((dipole.cross(z_axis).normalize()).dot(x_axis));
	      }
	      P1.add(dipole.dot(z_axis));
	      norm_P1.add((dipole.cross(z_axis).normalize()).dot(x_axis));
	    }
	    else {
	      in_out[i].push_back(false);
	    }
	      
	  }
	  output<< time;
	  output.precision(4);
	  output<<setw(15)<<Total_dipole.abs();
	  Total_dipole =Total_dipole.normalize();
	  dot_collective=Total_dipole.dot(z_axis);
	  P1_collective.add(dot_collective);
	  output<<setw(15)<<dot_collective;
	  Total_dipoles.push_back(Total_dipole);

	  norm_Total_dipole=(Total_dipole.cross(z_axis)).normalize();
          dot_collective=norm_Total_dipole.dot(x_axis);
          norm_P1_collective.add(dot_collective);
          output<<setw(15)<<dot_collective<<endl;
	  norm_Total_dipoles.push_back(norm_Total_dipole);
          loads_in_time.push_back(load_indexes);



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

    output.close();
    // full distribution single P1
    output.open("P1_full.dat");
    output<<"#Normalized P1 along tube of  : "<<atoms_count<<endl;
    output<<"#P1            Pdf"<<endl;
    P1.write_normalized(output);
    output.close();
   
    // full distribution single P1
    output.open("P1_full_norm.dat");
    output<<"#Normalized P1(x) along tube of  : "<<atoms_count<<endl;
    output<<"#P1            Pdf"<<endl;
    norm_P1.write_normalized(output);



    output.close();
    // full distribution collective P1
    output.open("P1_collective.dat");
    output<<"#Normalized collective P1 along tube of  : "<<atoms_count<<endl;
    output<<"#P1            Pdf"<<endl;
    P1_collective.write_normalized(output);
    output.close();
    


    output.close();
    // full distribution collective P1
    output.open("P1_collective_norm.dat");
    output<<"#Normalized collective P1 along tube of  : "<<atoms_count<<endl;
    output<<"#P1(x)            Pdf"<<endl;
    norm_P1_collective.write_normalized(output);
    output.close();

    
    // Distributions at each bin  of the pore and average P1 of each bin
    output.open("avg_P1_along_pore.dat");
    output<<"# <P1> parameter along pore axis of  : "<<atoms_count<<endl;
    output<<"# z         <P1>      error"<<endl;
    output.precision(4);
    for(int i = 0; i < hists1D.size(); i++) {
      P1.clear();
      Stat<double> P1_error;
      double z =(i)*bin_size + min_z +  bin_size*0.5;
      for(int j = 0; j < hists1D[i].size(); j++) {
	P1.add(hists1D[i][j]);
	P1_error.addval(hists1D[i][j]);
      }
      if (hists1D[i].size() > 0)
	output<<setw(15)<<z<<setw(15)<<P1_error.ave()<<setw(15)<<P1_error.ee()<<endl;
      // write P1 distributions for each bin
      if (WriteDistBin) {
	ofstream output2;
	string name1 = "z_";
	string name2 = static_cast<ostringstream*>( &(ostringstream() << z) )->str(); //convert int into string
	string name3 = "P1_" + name1 + "." + name2 + ".dat";
	output2.open(name3.c_str());
	output2<<"#Normalized P1 in "<<z<<" of : "<<atoms_count<<endl;
	output2<<"#P1              Pdf"<<endl;
	P1.write_normalized(output2);
	output2.close();
      }
      
    }
    output.close();


    //1d vector to auto-correlation 
    vector<double> P1_list(window,0.0);
    vector<double> norm_P1_list(window,0.0);
    // 1d vector to count for averaging
    vector<int> average_counter(window,0);
    
#ifdef OMP
#pragma omp parallel for num_threads(threads)
#endif
    for(int it = 0; it < frames-window; it++) { // loop trough origins
      int window_counter = 0;
      for(int j = it; j < it+window; j++) { //loop through frames starting from it
	P1_list[window_counter]+=(Total_dipoles[j]).dot((Total_dipoles[it]));
	norm_P1_list[window_counter]+=(norm_Total_dipoles[j]).dot((norm_Total_dipoles[it]));
	average_counter[window_counter]++;
	window_counter++;
      }
    }

    //output collective Dipole time-correlation data
   
    output.open("Mu_tcorr_collective.dat");
    output<<"#<coll_mu(t)dot coll_mu(t0) of atoms "<<atoms_count<<endl;
    output<<"#Time[ps]          <mu*mu>               <norm_mu*norm_mu>"<<endl;
    output.precision(4);
    
    for(int i = 0; i < P1_list.size(); i++) {
      //cout<<average_counter[i]<<endl;
      output<<setw(15)<<i*time.dt()<<setw(15)<<(P1_list[i]/average_counter[i])
	    <<setw(15)<<(norm_P1_list[i]/average_counter[i])<<endl;
    }
    output.close();
    


    
    //1d vector to auto-correlation 
    fill(P1_list.begin(),P1_list.end(),0.0);
    fill(norm_P1_list.begin(),norm_P1_list.end(),0.0);
    // 1d vector to count for averaging
    fill(average_counter.begin(),average_counter.end(),0);    

if (!FullCorrelation)
#ifdef OMP
#pragma omp parallel  num_threads(threads)
#endif
  for(int it = 0; it < frames-window; it++) {
	int window_counter = 0;
	double square_d = 0;
	for(int j = it; j < it+window; j++) { //loop through frames starting from it
	  for (int h =0; h < loads_in_time[it].size() ; h++) {
	    int i=loads_in_time[it][h];
	    square_d = 0;
	    if (in_out[i][j]) {
	      P1_list[window_counter]+=(dipoles[i][j]).dot((dipoles[i][it]));
	      norm_P1_list[window_counter]+=(norm_dipoles[i][j]).dot((norm_dipoles[i][it]));
	      average_counter[window_counter]++;
	    }
	  }	
	  window_counter++;
	}
      }
    //Full correlation code of all  molecules

 if (FullCorrelation)
#ifdef OMP
#pragma omp parallel  num_threads(threads)
#endif
      for (int k =0; k < at.size(); k+=NatMol) { //loop trough atoms
	for(int it = 0; it < frames-window; it++) { // loop trough origins
	  int window_counter = 0;
	  for(int j = it; j < it+window; j++) { //loop through frames starting from it
	    P1_list[window_counter]+=(dipoles[k][j]).dot((dipoles[k][it]));
	    average_counter[window_counter]++;
	    window_counter++;
	  }
	}
      }
 

    //output Dipole time-correlation data
   
    output.open("Mu_tcorr.dat");
    output<<"#<mu(t)dot mu(t0) of atoms "<<atoms_count<<endl;
    output<<"#Time[ps]          <mu*mu>               <norm_mu*norm_mu>"<<endl;
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





