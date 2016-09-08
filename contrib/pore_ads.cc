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
/* Program tha computes 2D histograms along the plane orthogonal to the pore axis i.e (x-y plane). The pore needs to be aligned in z.
   Additionally it computes a kind of rdf with respect to the x-y center of the pore, where the computed distance is only in the orh-t
   gonal plane (x-y) integrating over the pore main axis length. 
   The total adsorbed molecules based on the x-y distance to the center of the pore are also computed. employing a simple cutoff.
*/ 


const double PI =  3.14159265359;

//function to check normalize options
void check_norm_argument(string normalization);

//class declaration to unbin 1D and 2D histograms
class Unbin  {
public:
  Unbin(  std::vector<vector<int> >   & histogram, double  min_D1, double  min_D2,double  bin_1D, double bin_2D, double  normalization)
    : hist2D(histogram),min_1D(min_D1),min_2D(min_D2),bin_size_1D(bin_1D),bin_size_2D(bin_2D),norm_factor(normalization)
  {
  
  }
  
  Unbin(  std::vector<int>  &  histogram, double  min_D1, double  bin_1D, double  normalization)
    : hist1D(histogram),min_1D(min_D1),bin_size_1D(bin_1D),norm_factor(normalization)
  {
  
  }

  Unbin(  std::vector<int>  &  histogram, double  min_D1, double  bin_1D)
    : hist1D(histogram),min_1D(min_D1),bin_size_1D(bin_1D)
  {
  
  }



  ~Unbin() {}

  void unbin2D(std::ofstream & output);
  void unbin1D(std::ofstream & output);
  void unbin1D_rad (std::ofstream & output);
  void setnorm (double norm);
  
 private:
  std::vector<vector<int> > hist2D;
  std::vector<int> hist1D;
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
         <<"bins"<< "traj" << "ref" << "offset" << "cutoff" 
	 << "pore_radius"<<"normfact"<< "threads"<< "first"
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
  usage += "\t@cutoff   <adsorbtion max cutoff >\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@normfact   <normalization factor >\n";
  usage += "\t@threads   <number of threads  >\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    // get the @offset argument
    int bins = args.getValue<int>("bins");
    // get normalization argument
    string normalization;    
    {
      Arguments::const_iterator iter = args.lower_bound("normfact");
      Arguments::const_iterator to = args.upper_bound("normfact");
      for (; iter != to; iter++) {
	normalization = iter->second.c_str();
      }
    }
    
    
    check_norm_argument(normalization);
    

    // get the @offset argument
    double offset=0.0;
    if ( args.count("offset") > 0)
      offset=args.getValue<double>("offset");
    
    // get the @cutoff argument
    double cutoff = args.getValue<double>("cutoff");

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
    
    // 1d vector to store histogram in radial distances
    vector<int> hist1D(bins);  
    // 2d vector to store 2d histogram in x-y plane
    vector<vector<int> >  hist2D(bins,vector<int>(bins,0));
    //define rest of the variables
    double max_x,min_x, max_y,min_y,max_z,min_z,bin_size_x,bin_size_y,bin_size,center_x,center_y,rad2,rad;
    double step_x,step_y,step,rad_dist;
    Vec x_y(0.0, 0.0, 0.0);
    int x, y, bin;
    Vec vec_box;
    // total adsorbed molecules
    int total_ads = 0;
    ofstream output;
    output.open("adsorbed_TS.out");
    output << "# Total Pore adsorbed atoms of "<<atoms_count<<endl;
    output << "# Time[ps]   [N]"<<endl;
    
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
	  total_ads=0;
	  // do the fitting
	  (*pbc.*gathmethod)();
	  rf.fit(&sys);
	  //loop over pore molecules and get max and min z
	
	  min_z = pore_atoms.pos(0)[2];
	  max_z = pore_atoms.pos(0)[2];
	  		
	  vec_box    = (sys.box().K() + sys.box().L() + sys.box().M());
	
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
			
	  // set min and max for binnig
	  min_x = center_x - vec_box[0]*0.5;
	  max_x = center_x + vec_box[0]*0.5;
	  min_y = center_y - vec_box[1]*0.5;
	  max_y = center_y + vec_box[1]*0.5;
	
	  bin_size_x =(max_x-min_x)/(bins);
	  bin_size_y =(max_y-min_y)/(bins);
	
	  // set maximal radial distance in x-y plane
	  rad2= vec_box[0]*0.5*vec_box[0]*0.5 + vec_box[1]*0.5*vec_box[1]*0.5;
	  rad = sqrt(rad2);
	  bin_size = (rad)/(bins);
	
	  // loop through selected atoms

#ifdef OMP
#pragma omp parallel for  num_threads(threads)	\
  reduction(+: total_ads) private (x_y,rad_dist,step,bin,step_x,step_y,x,y)
#endif
	  for(int i = 0; i < count_atoms.size(); i++) {
	  
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
		step= floor( (rad_dist/bin_size));
		// Arrays elements start from zero
		bin=int(step);
		if (bin >= 0 && bin < bins)
		  hist1D[bin]++;
		//binnig 2d hist
		step_x= floor( (x_y[0]-min_x)/(bin_size_x));
		step_y= floor( (x_y[1]-min_y)/(bin_size_y));
		// Arrays elements start from zero
		x=int(step_x);
		y=int(step_y);
		// Reject anystep that is bigger than the bin size (avoids a seg fault) 
		if (x>=0 && x<bins && y>0 && y<bins)
		  hist2D[x][y]++;
		  
	      }		  
	    }
	  }
	  output << time;
	  output<<setw(9)<<total_ads<<endl;
	}
      }
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    output.close();
    double norm=1;
    double box_volume=1;
    double slice_volume=1;
    double Total_density,const_factor;
    
    if (normalization.compare("none") == 0) 
      norm=1;
    if (normalization.compare("average") == 0)
      norm=frames;
    if (normalization.compare("avgprob") == 0)
      norm=count_atoms.size()*frames;
    
    if (normalization.compare("rdf") == 0) {
      box_volume   = (vec_box[0])*(vec_box[1])*(vec_box[2]);
      slice_volume = bin_size_x*bin_size_y*(max_z-min_z);
      norm =(count_atoms.size()*frames*slice_volume)/box_volume;
    }
    
    // 2D histogram
        
    output.open("2Ddensity.out");
    output << "# 2D density across pore of "<<atoms_count<<endl;
    output << "# X   Y     [N]"<<endl;
    output.precision(3);
    
    Unbin Unbin_object_2D(hist2D,min_x,min_y,bin_size_x,bin_size_y,norm);
    Unbin_object_2D.unbin2D(output);
    output.close();
    
    
    // 1D histogram
    
    output.open("Radial_density.out");
    output << "# 1D (x-y plane) density across pore of "<<atoms_count<<endl;
    output << "# R    [N]"<<endl;   
    output.precision(3);

    Unbin Unbin_object_1D(hist1D,0.0,bin_size);
    if (normalization.compare("rdf") != 0) {
      Unbin_object_1D.setnorm(norm);
      Unbin_object_1D.unbin1D(output);
    }     
    if (normalization.compare("rdf") == 0) {
      Total_density = count_atoms.size()/box_volume;
      // const factor of the volume of a cylindrical shell 
      const_factor  = PI*(max_z-min_z)*bin_size*bin_size;
      // norm factor (dNs/Vshell)/(Ntotal/Vtotal)
      norm = Total_density*const_factor;
      // Include the average
      norm*=frames;
      Unbin_object_1D.setnorm(norm);
      Unbin_object_1D.unbin1D_rad(output);      
    } 
	
    output.close();
    
       
        
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

//functions and objects implementations

void check_norm_argument(string normalization) {
  string usage; 
  //  
  if ( normalization.compare("none") != 0  && normalization.compare("average") != 0 &&  normalization.compare("avgprob") != 0 && normalization.compare("rdf") != 0 ) {
    usage+= normalization;
    usage+= " option not known\n";
    usage+= "@normfact options are:\n";
    usage += "\tnone    = raw distribution in x-y plane/radial\n";
    usage += "\taverage = Time average of raw distribution in x-y plane/radial\n";
    usage += "\tavgprob = Time average of Nparticles(bin)/NtotParticles in x-y plane/radial\n";
    usage += "\trdf     = rdf type of normalization\n";
    cout<<usage;
    std::exit(EXIT_FAILURE); 
  } else
    cout<<"\n#Employing "<<normalization<<" normalization"<<endl;  
    
}


// Unbin class implementation
void Unbin::unbin2D(std::ofstream & output){
  double X,Y;
  for(int i = 0; i < hist2D.size(); i++) {
    X =(i)*bin_size_1D + min_1D + bin_size_1D*0.5;
    for(int j = 0; j < hist2D[0].size(); j++) {
      Y =(j)*bin_size_2D +min_2D + bin_size_2D*0.5;
      //cout<<hist2D[i][j]<<" "<<norm_factor<<endl;
      output<<setw(15)<<X<<setw(15)<<Y<<setw(15)<<(hist2D[i][j])/norm_factor<<endl;
    }
  }
}

void Unbin::unbin1D(std::ofstream & output){
  double X;
  for(int i = 0; i < hist1D.size(); i++) {
    X =(i)*bin_size_1D + min_1D + bin_size_1D*0.5;
    output<<setw(15)<<X<<setw(15)<<(hist1D[i])/norm_factor<<endl;
  }
}

// Unbins a cylindrically Radially binned histogram 
void Unbin::unbin1D_rad (std::ofstream & output) {
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

void Unbin::setnorm(double norm) {
  norm_factor =norm;
}


