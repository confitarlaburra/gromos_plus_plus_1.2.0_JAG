/**
 * @file pore_hsolv
 * Calculates the solvation energy of loaded particles
 * in a carbon nanotube
 */
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
#include <iomanip>  


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gmath/Physics.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/groTime.h"
// added by jag
#include "../src/gmath/Stat.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/args/GatherParser.h"


using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;
using namespace fit;

// Functions declarations

/* LJ = Computes lenard-jones interaction
   iac1 = atom 1
   iac2 = atom 2
   r2   = distance between atom 1 and atom 2
   gff  = force field
   Returns LJ interaction of  a given pair
*/
double LJ(int iac1, int iac2, double r2, GromosForceField &gff);


/* Coulomb = Computes Coulomb interaction
   c1 = atom 1
   c2 = atom 2
   r2 = distance between atom 1 and atom 2
   eps = permitivity of space (electric constant)
   kap = debye length
   cut = cutoff
   Returns reaction field approx. for the Coulombic interaction of a given pair 
*/
double Coulomb(double c1, double c2, double r2, double eps, double kap, double cut);


/* find_center = Finds the closest molecule to the center of the pore in the z axis 
   and  both particles at the pore mouths. 
   input:
   indexes array = array with molecules at pore at each frame.
   atoms = selected atoms
   Returns the index of atom closest to the center of the pore:
   access element: load_indexes[find_center(load_indexes,at,min_z,max_z)];
*/
int find_center( vector<int> &indexes, AtomSpecifier &atoms, double &min_z, double &max_z);


/* find_max = Finds the closest molecule to the upper mouth  of the pore in the z axis.
   input:
   indexes array = array with molecules at pore at each frame.
   atoms = selected atoms
   max_z = max coordinates
   Returns the index of atom closest to the center of the pore:
   access element: load_indexes[find_center(load_indexes,at,max_z)];
*/
int find_max( vector<int> &indexes, AtomSpecifier &atoms,double &max_z);


/* find_min = Finds the closest molecule to the upper mouth  of the pore in the z axis.
   input:
   indexes array = array with molecules at pore at each frame.
   atoms = selected atoms
   min_z = max coordinates
   Returns the index of atom closest to the center of the pore:
   access element: load_indexes[find_min(load_indexes,at,max_z)];
*/
int find_min( vector<int> &indexes, AtomSpecifier &atoms, double &min_z);

/*
  typedef to hold VdW and Coulomb time series
 */
typedef  map <string, Stat<double> > EnergyMap;

/*
  typedef to hold Components of the system
 */
typedef  map <string,EnergyMap > CompMap;



int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "rf" << "T" << "traj" << "time"
	 << "atoms"<<"ref"<<"offset"<<"pore_radius"<<"first"
	 << "last"<< "pore"<<"threads"<<"solvent"<<"NatMol"
	 <<"MaxLoad"<<"solv_atoms";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t[@pbc     <boundary type (read from GENBOX block if not specified)> [<gather method>]>]\n";
  usage += "\t@rf       <cut off radiuse> <epsilon> <kappa>\n";
  usage += "\t[@time    <time and dt, used to overwrite the time of read from the trajectory file>]\n";
  usage += "\t@T        <temperature, used for the RT volume expansion term>\n";
  usage += "\t@traj      <simulation trajectory or coordinate file>\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@NatMol   <Number of atoms of molecule>\n";
  usage += "\t@solvent   <Solvent atoms>\n";
  usage += "\t@solv_atoms   <Solvent is equivalent to atoms to follow (loaded)>\n 1:yes or 0:no";
  usage += "\t@last   <solvent atoms  >\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@offset   <offset for the main axis of the pore>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@threads    <thread number >";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  usage += "\t@MaxLoad   <Max load of pore with selected atoms>\n";
  

  try {
    Arguments args(argc, argv, knowns, usage);

    // handle the time
    Time time(args);
    
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
    
    int threads= 1;
    if ( args.count("threads") > 0)
      threads = args.getValue<int>("threads");
       
    // get NatmMol 
    int NatMol = args.getValue<int>("NatMol");

    // get NatmMol 
    int MaxLoad = args.getValue<int>("MaxLoad");

    // get NatmMol 
    bool solv_atoms = args.getValue<bool>("solv_atoms");
    
    // read the temperature
    if(args.count("T") != 1) {
      throw gromos::Exception(argv[0], "check @T argument to specify the temperature");
    }
    double T;
    {
      stringstream ss;
      ss << args["T"];
      ss >> T;
      if(ss.fail() || ss.bad()) {
        stringstream msg;
        msg << "could not convert " << args["T"] << " to be used as temperature";
        throw gromos::Exception(argv[0], msg.str());
      }
    }
    
    // read reaction field parameters
    double eps, kap, cut;
    if (args.count("rf") != 3) {
      throw gromos::Exception(argv[0], "too few arguments for the reaction-field parameters (@rf)");
    } else {
      Arguments::const_iterator it = args.lower_bound("rf");
      string ss;
      ss = it->second;
      cut  = atof(ss.c_str());
      ++it;
      ss = it->second;
      eps  = atof(ss.c_str());
      ++it;
      ss = it->second;
      kap  = atof(ss.c_str());
    }

   
    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System refSys(it.system());
    GromosForceField gff = it.forceField();
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
    
    
    //fit atoms (the pore is aligned in z)
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
      throw gromos::Exception("pore_energy",
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
      throw gromos::Exception("pore_energy",
			      "No atoms to calculate!");



     // get solvent atoms
    AtomSpecifier solvent(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("solvent");
      Arguments::const_iterator to = args.upper_bound("solvent");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        solvent.addSpecifier(spec);
      }
    }
    if ( !solvent.size() )
      throw gromos::Exception("pore_energy",
			      "No atoms to calculate pore dimensions!");
    
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);


    double  max_z,min_z,center_x,center_y;
    double rad2= pore_radius*pore_radius;
    /*
       multidimensional map structure eg:
       Ene3Dmap[A][B][C][D] is a stat vector of energies with:
       A = position un pore i.e. "mouth","center" and "full"
       B = load number i.e. 1 to n
       C = system component i.e. "Loaded", "Solvent"  "Pore" "Self" and "Total"
       D = nonbonded energy component i.e. "LJ" of "RF" and "Total"
    */
    map<string,vector<CompMap> > Ene3Dmap;
    //init multidimensional map
    {
      // init A  and B dimension of map
      Ene3Dmap["center"].resize(MaxLoad);
      Ene3Dmap["mouth"].resize(MaxLoad);
      Ene3Dmap["full"].resize(MaxLoad);
      // define a dummy vector to init C dimension
      vector<string> components(6);
      components[0]="Self";
      components[1]="Loaded";
      components[2]="Pore";
      components[3]="Solvent";
      components[4]="Notloaded";
      components[5]="Total";
      
      // define a dummy vector to init D dimension
      vector<string> ene_terms(3);
      ene_terms[0]="LJ";
      ene_terms[1]="RF";
      ene_terms[2]="Total";
      //define a dummy Stat<double> to init final dimension of map
      Stat<double> dummy_stat;
      //dummy_stat.addval(5);
      // deinfe iterator for init of mulidim map
      map<string,vector<CompMap> >::iterator Ene3Dmap_it;
      vector<string>::iterator ene_it;
      vector<string>::iterator comp_it;
      for (Ene3Dmap_it = Ene3Dmap.begin();
	   Ene3Dmap_it!= Ene3Dmap.end(); Ene3Dmap_it++) { // A dimension
	for (int i=0; i<Ene3Dmap[Ene3Dmap_it->first].size(); i++) { // B dimension
	  for (comp_it = components.begin();
	       comp_it != components.end(); comp_it++) { // C dimension
	    for (ene_it = ene_terms.begin();
		 ene_it != ene_terms.end(); ene_it++) { // D dimension
	      (((Ene3Dmap[Ene3Dmap_it->first])[i])[*comp_it])[*ene_it]  = dummy_stat;
	    }
	  }
	}
      }
	
    }
    


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
	  // save non fitted coordinates
	  // vectors to store non fitted postions
	  vector<Vec> at_NF (at.size());
	  vector<Vec> solvent_NF (solvent.size()); 
	  vector<Vec> pore_atoms_NF (pore_atoms.size());
	  
	  for(int i = 0; i < pore_atoms.size(); i++) {
	    pore_atoms_NF[i]=pore_atoms.pos(i);
	  }
	  for(int i = 0; i < solvent.size(); i++) {
	    solvent_NF[i]=solvent.pos(i);
	  }
	  for(int i = 0; i < at.size(); i++) {
	    at_NF[i]=at.pos(i);
	  }
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
	  vector<int> load_indexes,unload_indexes;
	  int load;
	  map <string,vector<int> > Tube_comp;
	  vector<int> dummy;
	  Tube_comp["center"]= dummy;
	  Tube_comp["full"]  = dummy;
	  Tube_comp["mouth"] = dummy;
	  int center,min,max,index;

	  // loop only considering center atoms
	  
	  for(int i = 0; i < at.size(); i+=NatMol) {
	    // if inside pore
	    if (((at.pos(i)[2] >= (min_z - offset)) && (at.pos(i)[2] <= (max_z + offset) ))  &&  
		( (at.pos(i)[0]-center_x)*(at.pos(i)[0]-center_x) + (at.pos(i)[1]-center_y)*(at.pos(i)[1]-center_y) <=  rad2) ) {
	      load_indexes.push_back(i); // only add center atom to count load
	      for (int j =0; j<NatMol; j++)
		Tube_comp["full"].push_back(i+j);
	    } else 
	      for (int j =0; j<NatMol; j++) 
		unload_indexes.push_back(i+j);
	  }
	  load = load_indexes.size();
	  if ( load > 0) {
	    center= load_indexes[find_center(load_indexes,at,min_z,max_z)];
	    min = load_indexes[find_min(load_indexes,at,min_z)];
	    max = load_indexes[find_max(load_indexes,at,max_z)];
	    for (int j =0; j<NatMol; j++) {
	      Tube_comp["center"].push_back(center+j);
	      Tube_comp["mouth"].push_back(min+j);
	    }
	    for (int j =0; j<NatMol; j++) {
	      Tube_comp["mouth"].push_back(max+j);
	    }
	  }
	  double RF_self,VdW_self,RF_loaded,VdW_loaded,RF_pore,VdW_pore,RF_solvent,VdW_solvent,RF_unload,VdW_unload,RF_total,VdW_total;
	  map<string,vector<int> >::iterator Tube_comp_it;
	  if ( load > 0 && load <= MaxLoad)
	    for (Tube_comp_it = Tube_comp.begin();
		 Tube_comp_it!= Tube_comp.end(); Tube_comp_it++) {
	      RF_self=0.0;
	      VdW_self=0.0;
	      RF_loaded=0.0;
	      VdW_loaded=0.0;
	      RF_pore=0.0;
	      VdW_pore=0.0;
	      RF_unload=0.0;
	      VdW_unload=0.0;
	      RF_solvent=0.0;
	      VdW_solvent=0.0;
	      RF_total=0.0;
	      VdW_total=0.0;
#ifdef OMP
#pragma omp parallel for reduction(+:VdW_self,RF_self,RF_loaded,VdW_loaded,RF_pore,VdW_pore,RF_unload,VdW_unload,RF_solvent,VdW_solvent) \
schedule(dynamic) num_threads(threads) 
#endif
	      for (int i=0; i<Tube_comp[Tube_comp_it->first].size(); i++) {
		Vec pos1 = at_NF[Tube_comp[Tube_comp_it->first][i]];
		int iac1       = at.iac(Tube_comp[Tube_comp_it->first][i]);
		double charge1 =at.charge(Tube_comp[Tube_comp_it->first][i]);
		int molNum1    = at.mol(Tube_comp[Tube_comp_it->first][i]);
		//self loop
		for (int j=i+1; j<Tube_comp[Tube_comp_it->first].size(); j++) {
		  int molNum2 = at.mol(Tube_comp[Tube_comp_it->first][j]);
		  if (molNum1  != molNum2 ) {
		    Vec pos2=at_NF[Tube_comp[Tube_comp_it->first][j]];
		    double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
		    int iac2=at.iac(Tube_comp[Tube_comp_it->first][j]);
		    double charge2 =at.charge(Tube_comp[Tube_comp_it->first][j]);
		    VdW_self += LJ(iac1, iac2, r2, gff);
		    RF_self  += Coulomb(charge1, charge2, r2, eps, kap, cut);
		  }
		}
		// loaded loop
		for (int j=0; j<Tube_comp["full"].size(); j++) {
		  int molNum2 = at.mol(Tube_comp["full"][j]);
		  if (molNum1  != molNum2 ) {
		    Vec pos2=at_NF[Tube_comp["full"][j]];
		    double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
		    int iac2=at.iac(Tube_comp["full"][j]);
		    double charge2 =at.charge(Tube_comp["full"][j]);
		    VdW_loaded += LJ(iac1, iac2, r2, gff);
		    RF_loaded  += Coulomb(charge1, charge2, r2, eps, kap, cut);
		  }
		}
		// Pore loop
		for (int j=0; j<pore_atoms.size(); j++) {
		  Vec pos2=pore_atoms_NF[j];
		  double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
		  int iac2=pore_atoms.iac(j);
		  double charge2 =pore_atoms.charge(j);
		  VdW_pore += LJ(iac1, iac2, r2, gff);
		  RF_pore  += Coulomb(charge1, charge2, r2, eps, kap, cut);
		}
		// solvent loop
		if(!solv_atoms)
		  for (int j=0; j<solvent.size(); j++) {
		    Vec pos2=solvent_NF[j];
		    double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
		    int iac2=solvent.iac(j);
		    double charge2 =solvent.charge(j);
		    VdW_solvent += LJ(iac1, iac2, r2, gff);
		    RF_solvent  += Coulomb(charge1, charge2, r2, eps, kap, cut);
		  }
		// unload loop
		for (int j=0; j<unload_indexes.size(); j++) {
		  Vec pos2=at_NF[unload_indexes[j]];
		  double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
		  int iac2=at.iac(unload_indexes[j]);
		  double charge2 =at.charge(unload_indexes[j]);
		  VdW_unload += LJ(iac1, iac2, r2, gff);
		  RF_unload  += Coulomb(charge1, charge2, r2, eps, kap, cut);
		}
	      }
	      // Accumulate energies
	      
	      //Self interactions
	      Ene3Dmap[Tube_comp_it->first][load-1]["Self"]["RF"].addval(RF_self);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Self"]["LJ"].addval(VdW_self);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Self"]["Total"].addval(RF_self+VdW_self);
	      RF_total += RF_self;
	      VdW_total+= VdW_self;
	      //loaded interactions
	      if (Tube_comp_it->first.compare("full") != 0) {
		Ene3Dmap[Tube_comp_it->first][load-1]["Loaded"]["RF"].addval(RF_loaded);
		Ene3Dmap[Tube_comp_it->first][load-1]["Loaded"]["LJ"].addval(VdW_loaded);
		Ene3Dmap[Tube_comp_it->first][load-1]["Loaded"]["Total"].addval(RF_loaded+VdW_loaded);
		RF_total += RF_loaded;
		VdW_total+= VdW_loaded;
	      }
	      else {
		Ene3Dmap[Tube_comp_it->first][load-1]["Loaded"]["RF"].addval(RF_loaded*0.5);
		Ene3Dmap[Tube_comp_it->first][load-1]["Loaded"]["LJ"].addval(VdW_loaded*0.5);
		Ene3Dmap[Tube_comp_it->first][load-1]["Loaded"]["Total"].addval((RF_loaded+VdW_loaded)*0.5);
	      }
	      //Pore interactions
	      Ene3Dmap[Tube_comp_it->first][load-1]["Pore"]["RF"].addval(RF_pore);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Pore"]["LJ"].addval(VdW_pore);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Pore"]["Total"].addval(RF_pore+VdW_pore);
	      RF_total += RF_pore;
	      VdW_total+= VdW_pore;
	      //Unloaded interactions
	      Ene3Dmap[Tube_comp_it->first][load-1]["Notloaded"]["RF"].addval(RF_unload);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Notloaded"]["LJ"].addval(VdW_unload);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Notloaded"]["Total"].addval(RF_unload+VdW_unload);
	      RF_total += RF_unload;
	      VdW_total+= VdW_unload;
	      //Solvent interactions= to unloaded interactions if solven is equal to loaded
	      if(!solv_atoms) {
		Ene3Dmap[Tube_comp_it->first][load-1]["Solvent"]["RF"].addval(RF_solvent);
		Ene3Dmap[Tube_comp_it->first][load-1]["Solvent"]["LJ"].addval(VdW_solvent);
		Ene3Dmap[Tube_comp_it->first][load-1]["Solvent"]["Total"].addval(RF_solvent+VdW_solvent);
		RF_total += RF_solvent;
		VdW_total+= VdW_solvent;
	      } else {
		Ene3Dmap[Tube_comp_it->first][load-1]["Solvent"]["RF"].addval(RF_unload);
		Ene3Dmap[Tube_comp_it->first][load-1]["Solvent"]["LJ"].addval(VdW_unload);
		Ene3Dmap[Tube_comp_it->first][load-1]["Solvent"]["Total"].addval(RF_unload+VdW_unload);
	      }
	      // Add totals
	      Ene3Dmap[Tube_comp_it->first][load-1]["Total"]["RF"].addval(RF_total);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Total"]["LJ"].addval(VdW_total);
	      Ene3Dmap[Tube_comp_it->first][load-1]["Total"]["Total"].addval(RF_total+VdW_total);
	      
	    
	    }
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    // Final ouput of Data:
    map<string,vector<CompMap> >::iterator Ene3Dmap_it;
    for (Ene3Dmap_it = Ene3Dmap.begin();
	 Ene3Dmap_it!= Ene3Dmap.end(); Ene3Dmap_it++) {
      string name = Ene3Dmap_it->first +".dat"; 
      ofstream output(name.c_str());
      output<<"#Average Non-bonded Energy for atoms "<<atoms_count<<" at "<<Ene3Dmap_it->first<<" of pore"<<endl;
      for (int i=0; i<Ene3Dmap[Ene3Dmap_it->first].size(); i++) {
	CompMap::iterator componentA;
	EnergyMap::iterator energyB;
	output<<"\n\n#Load "<<i+1<<endl;
	output << setw(16) <<"Component"<<setw(18)<<"<LJ>"<<setw(16)<<"error "<<setw(16)
	       <<"<RF>"<<setw(16)<<"error"<<setw(16)<<"<Total>"<<setw(16)<<"error"<<endl;
	for (componentA = Ene3Dmap[Ene3Dmap_it->first][i].begin();
	     componentA != Ene3Dmap[Ene3Dmap_it->first][i].end(); componentA++) {
	  output<<"\n"<< setw(16) << componentA->first;
	  for (energyB = Ene3Dmap[Ene3Dmap_it->first][i][componentA->first].begin();
	       energyB != Ene3Dmap[Ene3Dmap_it->first][i][componentA->first].end();
	       energyB++)    {
	    output<<setw(16)<<energyB->second.ave()<<setw(16)<<energyB->second.ee();
	  }
	}
      }
      output.close();
    }
  }
  catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

//END//
// functions implementations//

// Lennard Jones
double LJ(int iac1, int iac2, double r2, GromosForceField &gff) {
  double r6 = r2 * r2 * r2;
  AtomPair ap(iac1, iac2);
  double c12 = gff.ljType(ap).c12();
  double c6 = gff.ljType(ap).c6();
  return ( (c12 /r6) - c6 ) / r6;
}
// Reaction field aproximation for Coulombic interactions
double Coulomb(double c1, double c2, double r2, double eps, double kap, double cut) {
  if (c1 == 0 || c2 == 0) {
    return 0.0;
  } else {
    double Crf = ((2 - 2 * eps)*(1 + kap * cut) - eps * (kap * cut)*(kap * cut)) /
            ((1 + 2 * eps)*(1 + kap * cut) + eps * (kap * cut)* (kap * cut));
    
    double r = std::sqrt(r2);
    double cut3 = cut * cut * cut;
    //return Crf;
    return (c1 * c2 * physConst.get_four_pi_eps_i() * ((1 / r) - (0.5 * Crf * r2 / cut3) - (1 - 0.5 * Crf / cut)));
  }
}

// Accesory functions
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

int find_max( vector<int> &indexes, AtomSpecifier &atoms, double &max_z) {
  double z  = 0.0;
  double small = abs ( atoms.pos(indexes[0])[2]-max_z);
  int index = 0;
  for (int i = 0; i < indexes.size(); i++) {
    z=abs(atoms.pos(indexes[i])[2]-max_z);
    if (z<small) {
      small=z;
      index =i;
    }
  }
  return index;
}

int find_min( vector<int> &indexes, AtomSpecifier &atoms, double &min_z) {
  double z  = 0.0;
  double small = abs ( atoms.pos(indexes[0])[2]-min_z);
  int index = 0;
  for (int i = 0; i < indexes.size(); i++) {
    z=abs(atoms.pos(indexes[i])[2]-min_z);
    if (z<small) {
      small=z;
      index =i;
    }
  }
  return index;
}
