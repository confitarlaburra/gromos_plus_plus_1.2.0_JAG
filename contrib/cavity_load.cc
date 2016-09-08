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
#include "delaunay.cc"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;


//JAG 2.6.2013
/* Cavity load counter using delaunay triangulation
*/ 


int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "cavity" << "atoms"  
         << "traj" << "outdelaunay";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@cavity    <cavity atoms>\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@outdelaunay   <out delaunay triangles in pdb format>\n";
  

  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    
    //read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    InG96 ic;
    if (args.count("traj") > 0) {
      ic.open(args.lower_bound("traj")->second);
      ic.select("ALL");
      ic >> refSys;
      ic.close();
    } else {
      throw gromos::Exception("frameout", "no trajectory specified (@traj)");
    }

   
    // The current system
    System sys(refSys);
    
    // get pore atoms
    AtomSpecifier cavity_atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("cavity");
      Arguments::const_iterator to = args.upper_bound("cavity");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        cavity_atoms.addSpecifier(spec);
      }
    }
    if ( !cavity_atoms.size() )
      throw gromos::Exception("cavity_load",
			      "No atoms to calculate cavity dimensions!");
    
        
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
      throw gromos::Exception("cavity_load",
			      "No atoms to calculate cavity loadings!");

    // Radius for spherical projection
    double R = 0.5;
    //matrix that contains shared triangles (init each element to 0)
    vector<vector<int> > shared;
    std::vector<int> Row(cavity_atoms.size(),0);
    for(int i = 0; i < cavity_atoms.size(); i++)
      shared.push_back(Row);
    // container of delaunay triangles
    vector <Delaunay> Delaunays;
    vector<Vec> old_positions;
    
    // output streams
    ofstream cl;
    cl.open("cavity_load.out");
    cl << "# Time series of cavity loadings"<<endl;
    cl << "# Time           Load"<<endl;
    ofstream out_pdb;
    if (args.count("outdelaunay") > 0)
      out_pdb.open("delaunays.pdb");
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);
       
    // loop over all trajectories 
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    for (; iter != to; iter++) {
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()) {
	ic >> sys >> time;
	// Shift to cog of cavity atoms
        PositionUtils::shiftToCog(&sys,cavity_atoms);
	//Sphere projection
	for(int i = 0; i < cavity_atoms.size(); i++) {
	  old_positions[i]=cavity_atoms.pos(i);
	  cavity_atoms.pos(i) *= R/(cavity_atoms.pos(i).abs2());
	}
	// Loop trough all cavity atoms and get all triplets
	Vec COG = PositionUtils::cog(sys,cavity_atoms); //should be zero
	for(int i = 0; i < cavity_atoms.size(); i++) {
	  for(int j = i+1; j < cavity_atoms.size(); j++) {
	    for(int k = i+2; k < cavity_atoms.size(); k++) {
	      if(shared[k][j]!=2) {
		Delaunay triangle(&cavity_atoms,i,j,k);
		triangle.set_cog();
		triangle.set_angle(pbc,sys,COG);
		//Delaunays.push_back(triangle);
		shared[i][j]++;
		shared[i][k]++;
		shared[j][k]++;
		// Loop trough all atoms and get delaunay triangles
		for(int l = 0; l < cavity_atoms.size(); l++) {
		  if (triangle.Calc_Angle(pbc,sys, COG, l) >= triangle.get_angle()) {
		    Delaunays.push_back(triangle);
		    break;
		  }
		}
	      }	    
	    }
	  } 
	}
	// project back to original shape
	for(int i = 0; i < cavity_atoms.size(); i++) {
	  cavity_atoms.pos(i)= old_positions[i] ;
	} 
	// print output pdb
	if (args.count("outdelaunay") > 0) {
	  out_pdb << "REMARK " << " Time" << time.time()<<endl;
	  for(int i = 0; i < Delaunays.size(); i++) {
	    Delaunays[i].writepdb(out_pdb,i) ;
	    out_pdb <<"ENDMDL"<<endl;
	  }
	}
	// 3.6.13 (check if actually works)


      }
    }
    cl.close();
    out_pdb.close();  
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
    return 0;
}
