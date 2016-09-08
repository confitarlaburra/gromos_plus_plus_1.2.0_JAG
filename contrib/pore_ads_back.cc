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


//JAG 23.4.2013
/* Simple porgram that count molecules within a rigid pore, like a carbon nanotube or an Aquaporin.
   Currently, it requires the reference frame to have the main axis of the pore aligned in the z-axis
*/ 



//template class to unbin 1D and 2D histograms
template <typename T1, typename T2> class Unbin  {
public:
  Unbin( const std::vector<vector<int> > & histogram, T1  min_D1, T1  min_D2,T1  bin_1D, T1  bin_2D, T2  normalization){
    set_hist (histogram);
    set_min_1D (min_D1);
    set_min_2D (min_D2);
    set_bin_size_1D(bin_1D);
    set_bin_size_2D(bin_2D);
    set_norm_factor(normalization);
  }


  Unbin( const std::vector<int> & histogram, T1  min_D1, T1  bin_1D, T2  normalization){
    set_hist (histogram);
    set_min_1D (min_D1);
    set_bin_size_1D(bin_1D);
    set_norm_factor(normalization);
  }

  Unbin();

  //~Unbin();

  void unbin2D(std::ofstream & output){
    for(int i = 0; i < hist2D.size(); i++) {
      double X =(i)*bin_size_1D + min_1D;
      for(int j = 0; j < hist2D[0].size(); j++) {
	double Y =(j)*bin_size_2D +min_2D;
	output<<setw(9)<<X<<setw(9)<<Y<<setw(9)<<(hist2D[i][j])/norm_factor<<endl;
      }
    }
  }

  void unbin1D(std::ofstream & output){
    for(int i = 0; i < hist1D.size(); i++) {
      double X =(i)*bin_size_1D + min_1D;
      output<<setw(9)<<X<<setw(9)<<(hist1D[i])/norm_factor<<endl;
    }
  }



  // 2d hist
  void set_hist (const std::vector<vector<int> > & histogram) {
    hist2D = histogram;
  }
  // 1D hist
  void set_hist (const std::vector<int> & histogram) {
    hist1D = histogram;
  }

  void set_min_1D (T1 min) {
    min_1D = min;
  }

  void set_min_2D (T1 min) {
    min_2D = min;
  }

  void set_bin_size_1D (T1 bin_1D) {
    bin_size_1D = bin_1D;
  }

  void set_bin_size_2D (T1 bin_2D) {
    bin_size_2D = bin_2D;
  }

  void set_norm_factor (T2 normalization) {
    norm_factor = normalization;
  }

private:
  std::vector<vector<int> > hist2D;
  std::vector<int> hist1D;
  T1 min_1D;
  T1 min_2D;
  T1 bin_size_1D;
  T1 bin_size_2D;
  T2 norm_factor;
  };
// end template class


int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "pore" << "atoms"  
         <<"bins"<< "traj" << "ref" << "offset" << "cutoff" 
	 << "pore_radius";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@bins  <Number of bins for the 2D histogram in X,Y>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@cutoff   <adsorbtion max cutoff >\n";
  usage += "\t@pore_radius   <pore radius >\n";
  

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
    
    // get the @cutoff argument
    double cutoff = args.getValue<double>("cutoff");

    // get the @pore_radius argument
    double pore_radius = args.getValue<double>("pore_radius");
    
      
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
      // Parse boundary conditions
      //Boundary *pbc = BoundaryParser::boundary(refSys, args);
      //parse gather method
      //Boundary::MemPtr gathmethod = args::GatherParser::parse(refSys, refSys, args);
      //(*pbc.*gathmethod)();
      
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
    
    // 1d vector to store histogram in radial distances
    vector<int> hist1D(bins+1,0);  
    // 2d vector to store 2d histogram in x-y plane
    vector<vector<int> >  hist2D(bins+1,vector<int>(bins+1,0));
    // Mass of selected particle
    double mass = count_atoms.mass(0);
    double max_x,min_x, max_y,min_y,bin_size_x,bin_size_y,bin_size,slice_volume;
    // total adsorbed molecules
    int total_ads = 0;
    ofstream adsorbed;
    adsorbed.open("adsrbed_TS.out");
    adsorbed << "# Total Pore adsorbed atoms of "<<atoms_count<<endl;
    adsorbed << "# Time[ps]   [N]"<<endl;

    
    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    int frames=0;
    
    //double bin_size;
    for (; iter != to; iter++) {
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()) {
	ic >> sys >> time;
	frames++;
	
	// do the fitting
	(*pbc.*gathmethod)();
	rf.fit(&sys);
	//loop over pore molecules and get max and min z
	
	//double min_z,max_z;
	double min_z = pore_atoms.pos(0)[2];
	double max_z = pore_atoms.pos(0)[2];
	min_x = count_atoms.pos(0)[0];
	max_x = count_atoms.pos(0)[0];
		
	double center_x,center_y;
	Vec vec_box    = (sys.box().K() + sys.box().L() + sys.box().M());
	
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
			
	// set min and max for binnig
	///*
	min_x = center_x - vec_box[0]*0.5;
	max_x = center_x + vec_box[0]*0.5;
	min_y = center_y - vec_box[1]*0.5;
	max_y = center_y + vec_box[1]*0.5;
	//*/
	bin_size_x =(max_x-min_x)/(bins);
	bin_size_y =(max_y-min_y)/(bins);
	slice_volume = bin_size_x*bin_size_y*(max_z-min_z);
	double rad2= vec_box[0]*0.5*vec_box[0]*0.5 + vec_box[1]*0.5*vec_box[1]*0.5;
	// max for radial distance binnig
	double rad = sqrt(rad2);
	bin_size = (rad)/(bins);
	
	// loop through selected atoms
		
	for(int i = 0; i < count_atoms.size(); i++) {
	  double step_x,step_y,step,rad_dist;
	  double new_x, new_y;
	  Vec x_y(0.0, 0.0, 0.0);
	  int X, Y, BIN;
	  
	  if ((count_atoms.pos(i)[2] >= (min_z - offset)) && (count_atoms.pos(i)[2] <= (max_z + offset) ) ) {
	    //scale coordinates relative to pore COG
	    x_y[0] = count_atoms.pos(i)[0] - center_x; 
	    x_y[1]=  count_atoms.pos(i)[1] - center_y;
	    // select molecules within a max radius 
	    if( ( x_y.abs2()) <=  rad2 ) {
	      rad_dist = x_y.abs();
	      if (rad_dist > pore_radius && rad_dist <= cutoff)
		total_ads++;
	      
	      //binnig 1D radial distance
	      step= floor( (rad_dist/bin_size) +1.0);
	      // Arrays elements start from zero
	      BIN=int(step)-1;
	      if (BIN>0 && BIN<bins)
		  hist1D[BIN]++;
	      //binnig 2d hist
	      step_x= floor( (x_y[0]-min_x)/(bin_size_x) +1.0 );
	      step_y= floor( (x_y[1]-min_y)/(bin_size_y) +1.0 );
	      // Arrays elements start from zero
	      X=int(step_x)-1;
	      Y=int(step_y)-1;
	      // Reject anystep that is bigger than the bin size (avoids a seg fault) 
	      if (X>0 && X<bins && Y>0 && Y<bins)
		hist2D[X][Y]++;
	    }
	  }
	}
      }
    }


    
    // 2D histogram
    ofstream density2D;
    density2D.open("2Ddensity.out");
    density2D << "# 2D density across pore of "<<atoms_count<<endl;
    density2D << "# X   Y     [N]"<<endl;
    
    density2D.precision(3);
    double norm =(frames*slice_volume)/mass;
    norm =1;
    Unbin<double,double> Unbin_object_2D(hist2D,min_x,min_y,bin_size_x,bin_size_y,norm);
    Unbin_object_2D.unbin2D(density2D);
    density2D.close();



    // 1D histogram
    ofstream density1D;
    density1D.open("Radial_density.out");
    density1D << "# 1D density across pore of "<<atoms_count<<endl;
    density1D << "# R    [N]"<<endl;   
    density1D.precision(3);
    //double norm =(frames*slice_volume)/mass;
    norm =1;
    double min=0.0;
    Unbin<double,double> Unbin_object_1D(hist1D,0.0,bin_size,norm);
    Unbin_object_1D.unbin1D(density1D);
    density1D.close();




    /*
    // Total Adsorbed time-series
    ofstream Adsorb_ts;
    Adsorb_ts.open("Adsorb_ts.out");
    Adsorb_ts << "# Time series of adsorbed "<<atoms_count<<endl;
    Adsorb_ts << "# Time   [N]"<<endl;
    */
     // output streams
   
    
        
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







