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


//JAG 23.4.2013
/* Simple porgrams that count molecules within a rigid pore, like a carbon nanotube or an Aquaporin.
   Currently, it requires the reference frame to have the main axis of the pore aligned in the z-axis
*/ 

// smoothing function that reduces fluctuations around pore boundaries
double smooth (double z, double r, double length, double z0, double r0);

int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "pore" << "atoms"  
         << "radius" << "traj" << "ref" << "offset" << "threads"
	 << "r0"<<"z0";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@radius  <pore radius>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@threads  <number of threads for openmp>\n";
  usage += "\t@r0  <sampled radius per loaded atom for smoothing>\n";
  usage += "\t@z0  <sampled axial  length perl loaded atom for smoothing>\n";
  

  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    // get the @radius argument
    double rad=args.getValue<double>("radius");
    double rad2=rad*rad;
    double center_x,center_y;
    
    int threads;
    if (args.count("threads") > 0) 
      threads=args.getValue<int>("threads");
    else
      threads=1;

    // get the @offset argument
    double offset=0.0;
    if ( args.count("offset") > 0)
      offset=args.getValue<double>("offset");
    double offset2=offset*offset;

    double z0=1;
    if ( args.count("z0") > 0)
      z0=args.getValue<double>("z0");
    
    double r0=1;
    if ( args.count("r0") > 0)
      r0=args.getValue<double>("r0");

        


      
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
      throw gromos::Exception("pore_loading",
			      "No atoms to calculate pore dimensions!");
    
        
    //get atoms to count
    AtomSpecifier count_atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        count_atoms.addSpecifier(spec);
      }
    }
    if ( !count_atoms.size() )
      throw gromos::Exception("pore_loading",
			      "No atoms to calculate pore loadings!");

    
    // output stream
    ofstream pp;
    pp.open("poreloadings.out");
    pp << "# Time series of Pore loadings with offset "<< offset <<endl;
    pp << "# Time           Load"<<endl;
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

       
    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    for (; iter != to; iter++) {
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()) {
	ic >> sys >> time;
	// do the fitting (it also center sytems on ref molecule)
	(*pbc.*gathmethod)();
	rf.fit(&sys);        
	//loop over pore molecules and get max and min z
	// init min_z and max_z to the max and min values of the box
       	double min_z = sys.box().Z();
	double max_z = -sys.box().Z();
	
	// Loop trough all pore atoms to get min and max z of the pore
	for(int i = 0; i < pore_atoms.size(); i++) {
	  center_x+=pore_atoms.pos(i)[0];
	  center_y+=pore_atoms.pos(i)[1];
	  if(pore_atoms.pos(i)[2] > max_z){
	    max_z = pore_atoms.pos(i)[2];
	  }
	  if(pore_atoms.pos(i)[2] < min_z){
	    min_z = pore_atoms.pos(i)[2];
	  }
	}
	
	center_x /=pore_atoms.size();
	center_y /=pore_atoms.size();
	
	double Pore_Length = max_z -min_z; 
		 
	// loop over atoms to count and check if they are within the pore
	double pore_wt_atoms=0.0;
	Vec x_y(0.0, 0.0, 0.0);
	Vec x_y_z(0.0, 0.0, 0.0);

#ifdef OMP
#pragma omp parallel for  num_threads(threads) \
  reduction(+: pore_wt_atoms) private(x_y,x_y_z)
#endif
	for(int i = 0; i < count_atoms.size(); i++) {
	  if ((count_atoms.pos(i)[2] >= (min_z)) && (count_atoms.pos(i)[2] <= (max_z) ) ) {
	    x_y[0]=count_atoms.pos(i)[0]-center_x; 
	    x_y[1]=count_atoms.pos(i)[1]-center_y;
	    if( ( x_y.abs2()) <=  rad2 )
	      pore_wt_atoms++;
	  }

	  //Apply smoothing function around pore mouths (offset radius)
	  
	  //Upper mouth
	  if  ( count_atoms.pos(i)[2] > max_z ) {
	    x_y_z[0]=count_atoms.pos(i)[0]-center_x; 
	    x_y_z[1]=count_atoms.pos(i)[1]-center_y;
	    x_y_z[2]=count_atoms.pos(i)[2]-max_z;
	    if(  x_y_z.abs2() <=  offset2 ) {
	      double z =count_atoms.pos(i)[2];
	      x_y_z[2]=0;// set it to zero so we only consider radial distance
	      pore_wt_atoms+=smooth (z, x_y_z.abs(),Pore_Length ,z0, r0);
	    }
	      
	  }

	  //Lower  mouth
	  if  ( count_atoms.pos(i)[2] < min_z ) {
	    x_y_z[0]=count_atoms.pos(i)[0]-center_x; 
	    x_y_z[1]=count_atoms.pos(i)[1]-center_y;
	    x_y_z[2]=count_atoms.pos(i)[2]-min_z;
	    if(  x_y_z.abs2() <=  offset2 ) {
	      double z =count_atoms.pos(i)[2];
	      x_y_z[2]=0;// set it to zero so we only consider radial distance
	      pore_wt_atoms+=smooth (z, x_y_z.abs(),Pore_Length ,z0, r0);
	    }
	      
	   }
	}
	// Round to the closest integer
	pore_wt_atoms=floor(pore_wt_atoms + 0.5);
	pp.precision(3);
	pp << time;
	pp.precision(3);
        pp << setw(15) <<  pore_wt_atoms << endl;
      }
    }
    pp.close();
  
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


// function declarations

double smooth (double z, double r,double length, double z0, double r0) {
  return exp( -(((abs(z)-length/2)/z0)*((abs(z)-length/2)/z0)) - (r/r0)*(r/r0) );
  
}
