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


//JAG 14.1.2014
/* Program tha computes the pair correlation function along the pore axis
*/ 


//class declaration to unbin 1D and 2D histograms
template <class T>
class Unbin  {
public:
  Unbin(  std::vector<vector<T> >   & histogram, double  min_D1, double  min_D2,double  bin_1D, double bin_2D, double  normalization)
    : hist2D(histogram),min_1D(min_D1),min_2D(min_D2),bin_size_1D(bin_1D),bin_size_2D(bin_2D),norm_factor(normalization)  {  
  }
  Unbin(  std::vector<vector<T> >   & histogram, double  min_D1, double  min_D2,double  bin_1D, double bin_2D)
    : hist2D(histogram),min_1D(min_D1),min_2D(min_D2),bin_size_1D(bin_1D),bin_size_2D(bin_2D){
  }
  Unbin(  std::vector<T>  &  histogram, double  min_D1, double  bin_1D, double  normalization)
    : hist1D(histogram),min_1D(min_D1),bin_size_1D(bin_1D),norm_factor(normalization){
  }
  Unbin(  std::vector<T>  &  histogram, double  min_D1, double  bin_1D)
    : hist1D(histogram),min_1D(min_D1),bin_size_1D(bin_1D){
  }
  ~Unbin() {}


  void unbin2D(std::ofstream & output);
  void unbin2D_rad(std::ofstream & output);
  void unbin1D(std::ofstream & output);
  void unbin1D_rad (std::ofstream & output);
  void setnorm (double norm);

 private:
  std::vector<vector<T> > hist2D;
  std::vector<T> hist1D;
  double min_1D;
  double min_2D;
  double bin_size_1D;
  double bin_size_2D;
  double norm_factor;
  };
// end  class


// main

int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "pore" << "atoms"  
         <<"bins"<< "traj" << "ref" << "offset" 
	 << "pore_radius"<< "threads"<< "first"
	 <<"last";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@bins  <Number of bins for each dimension in X,Y>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@threads   <number of threads  >\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    // get the @offset argument
    int bins = args.getValue<int>("bins");
    
    // get the @offset argument
    double offset=0.0;
    if ( args.count("offset") > 0)
      offset=args.getValue<double>("offset");
    
    
     // get the @first argument
    int first = args.getValue<int>("first");
    
    // get the @first argument
    int last = args.getValue<int>("last");

    // get the @pore_radius argument
    double pore_radius = args.getValue<double>("pore_radius");
    //get threads
    int threads;
    if (args.count("threads") > 0) 
      threads=args.getValue<int>("threads");
    else
      threads=1;      

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
    // 1d vector to store histogram  alonn pore axis
    vector<double> hist1D(bins);
    //define rest of the variables
    double max_x,min_x, max_y,min_y,max_z,min_z,bin_size_z,center_x,center_y,rad2,pore_length;
    double step_z;
    Vec x_y(0.0, 0.0, 0.0);
    int z;
    Vec vec_box;
    ofstream output;
        
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

	  pore_length = abs(max_z - min_z)+ 2*offset;
	 	
	  bin_size_z =pore_length/bins;
	  
	  // set maximal radial distance in x-y plane
	  rad2= pore_radius*pore_radius;
	  	
	  vector<int>load_indexes;

	  // loop through selected atoms

#ifdef OMP
#pragma omp parallel for  num_threads(threads)	\
  private (x_y)
#endif
	  for(int i = 0; i < count_atoms.size(); i++) {
	    if ((count_atoms.pos(i)[2] >= (min_z - offset)) && (count_atoms.pos(i)[2] <= (max_z + offset) ) ) {
	      //scale coordinates relative to pore COG
	      x_y[0] = count_atoms.pos(i)[0] - center_x; 
	      x_y[1]=  count_atoms.pos(i)[1] - center_y;
	      // select molecules within a max radius 
	      if( ( x_y.abs2()) <=  rad2 ) {
		load_indexes.push_back(i);
	      }
	    }
	  }

	  double inverse = (1/float(load_indexes.size())); 
#ifdef OMP
#pragma omp parallel for  num_threads(threads)	\
  private (step_z,z)
#endif

	  //loop trough loaded partilces and compute axial pair correlation function (bining)
	  
	  for(int i = 0; i < load_indexes.size(); i++) {
	    int atom_i=load_indexes[i];
	    //cout<<"i "<<count_atoms.pos(atom_i)[2]<<endl;
	    for(int j =0 ; j < load_indexes.size(); j++) {
	      if (j != i) {
		//binnig 1D axial distances
		int atom_j=load_indexes[j];
		//cout<<"j "<<count_atoms.pos(atom_j)[2]<<"delta "<<count_atoms.pos(atom_i)[2]-count_atoms.pos(atom_j)[2]<<endl;
		double distance = abs(count_atoms.pos(atom_i)[2]-count_atoms.pos(atom_j)[2]);
		step_z= floor(distance/(bin_size_z));
		//cout<<distance<<endl;
		z=int(step_z);
		// Arrays elements start from zero
		if (z >= 0 && z < bins)
		  hist1D[z]+=inverse ; //loads are not always the same thus we readily normalize by the load size
	      }
	    }
	  }		  
	

	}
      }
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    
    
    // 1D histogram
    
    output.open("Axial_pair_correlation.dat");
    output << "#axial distance across pore of "<<atoms_count<<endl;
    output << "# Axial rij    [N]"<<endl;   
    output.precision(3);

    Unbin<double> Unbin_object_1D(hist1D,0,bin_size_z);
    Unbin_object_1D.setnorm(frames);
    Unbin_object_1D.unbin1D(output);     
    	
    output.close();
    
       
        
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

//functions and objects implementations


template <class T>
void Unbin<T>::unbin2D(std::ofstream & output){
  double X,Y;
  for(int i = 0; i < hist2D.size(); i++) {
    X =(i)*bin_size_1D + min_1D;
    for(int j = 0; j < hist2D[0].size(); j++) {
      Y =(j)*bin_size_2D +min_2D + bin_size_2D*0.5 ;
      //cout<<hist2D[i][j]<<" "<<norm_factor<<endl;
      output<<setw(15)<<X<<setw(15)<<Y<<setw(15)<<(hist2D[i][j])/norm_factor<<endl;
    }
  }
}

// Cylindrical and axial histogram. Radial dimension is the second dimension of 2D array
template <class T>
void Unbin<T>::unbin2D_rad(std::ofstream & output){
  double X,Y,radial_norm_fact;
  int shell;
  for(int i = 0; i < hist2D.size(); i++) {
    X =(i)*bin_size_1D + min_1D;
    for(int j = 0; j < hist2D[0].size(); j++) {
      Y =(j)*bin_size_2D +min_2D + bin_size_2D*0.5;
      shell =2*(j+1) - 1;
      radial_norm_fact=norm_factor*shell*bin_size_1D;
      //cout<<hist2D[i][j]<<" "<<norm_factor<<endl;
      output<<setw(15)<<X<<setw(15)<<Y<<setw(15)<<(hist2D[i][j])/norm_factor<<endl;
    }
  }
}

template <class T>
void Unbin<T>::unbin1D(std::ofstream & output){
  double X;
  for(int i = 0; i < hist1D.size(); i++) {
    X =(i)*bin_size_1D + min_1D+ bin_size_1D*0.5;
    output<<setw(15)<<X<<setw(15)<<(hist1D[i])/norm_factor<<endl;
  }
}

template <class T> 
void Unbin<T>::unbin1D_rad (std::ofstream & output) {
  double radial_norm_fact,X;
  int shell;
  for(int i = 0; i < hist1D.size(); i++) {
    X =(i)*bin_size_1D + min_1D + bin_size_1D*0.5;
    //extra factor  for cylindrical shell i.e. each shell is larger in volume
    shell =2*(i+1) - 1;
    radial_norm_fact=norm_factor*shell;
    output<<setw(15)<<X<<setw(15)<<(hist1D[i])/radial_norm_fact<<endl;
  }
}
template <class T>
void Unbin<T>::setnorm(double norm) {
  norm_factor =norm;
}
