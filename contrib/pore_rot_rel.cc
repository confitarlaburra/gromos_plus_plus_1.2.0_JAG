/**
 * @file rot_rel.cc
 * Calculates autocorrelation functions for rotational relaxation times of particles loaded in zero centered cylindrical pore aligned in the z axis
 */

/**
 * The rotational relaxation time of molecules can be estimated from the
 * autocorrelation function of the Legendre polynomials of molecular axes 
 * @f$\vec{r}_i,\vec{r}_j$ and $\vec{k}@f$.
 * * @f[ C_1(t) = \left<\vec{r_i}(\tau) \cdot \vec{r_i}(\tau+t)\right>_{\tau} @f]
 * @f[ C_2(t) = \frac{1}{2} ( 3 \left<\vec{r_i}(\tau) \cdot \vec{r_i}(\tau+t)\right>^2 _{\tau} - 1 ) @f]
 *
 * Program rot_rel calculates the first and second order Legendre polynomials
 * and calculates the time correlation functions. The user specifies two of the
 * molecular axes, the third is defined as the cross product of the first two.
 * The program can average the correlation functions over multiple molecules in
 * the system using the flags @average
 * Note that the output of this program can also be produced by a combination
 * of programs @ref tser and @ref tcf.
 * This program is parallelised.
 *
 */

#ifdef OMP
#include <omp.h>
#endif
#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Value.h"
#include "../src/utils/VectorSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
//Added by JAG
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"


using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "ax1" << "ax2" << "traj" << "average" << "molecules"
	 <<"ref"<<"pore_radius"<<"pore"<<"first"<<"last"<<"window"<<"origin"<<"threads"<<"atoms"
	 <<"threads"<<"IN";

  string usage = argv[0];
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@ax1     <vector specifying molecular axis 1>\n";
  usage += "\t@ax2     <vector specifying molecular axis 2>\n";
  usage += "\t@average (average over all molecules)\n";
  usage += "\t@molecules <molecule numbers>(specify molecules for averaging)\n";
  usage += "\t@traj    <trajectory files>\n";
  usage += "\t@pore_radius   <pore radius >\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@window    <time window>\n";
  usage += "\t@threads    <thread number>\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@IN  <compute for inside pore (1) or outside pore >\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
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
	throw gromos::Exception("pore_rot_rel", "no trajectory specified (@traj)");
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
    string atoms_count;
    AtomSpecifier pore_atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("pore");
      Arguments::const_iterator to = args.upper_bound("pore");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atoms_count = spec;
	pore_atoms.addSpecifier(spec);
      }
    }
    if ( !pore_atoms.size() )
      throw gromos::Exception("pore_rot_rel",
			      "No atoms to calculate pore dimensions!");
    

    //get atoms to count
    
    AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
	
        at.addSpecifier(spec);
      }
    }
    if ( !at.size() )
      throw gromos::Exception("pore_rot_rel",
			      "No atoms to calculate!");
    // get simulation time
    Time time(args);
     
   // get the @first argument
    int first = args.getValue<int>("first");
    
    // get the @last argument
    int last = args.getValue<int>("last");
    
    // get the @pore_radius argument
    double pore_radius = args.getValue<double>("pore_radius");
   
     // get the @window argument
    int window = args.getValue<int>("window");

    // get the @threads argument
    int threads = args.getValue<int>("threads");

    // get the @threads argument
    bool IN = true;
    IN = args.getValue<bool>("IN");
    
        
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // read in the axes to follow
    // In this implementation, the user can specify two axes through a 
    // VectorSpecifier, which is either the vector connecting two atoms, or the
    // position of one atom.The third axis is defined as the cross-product of
    // these two.

    // vector of VectorSpecifier
    vector<utils::VectorSpecifier> vs1;
    vector<utils::VectorSpecifier> vs2;

    
    string s1=args["ax1"];
    string s2=args["ax2"];
    

    bool Va1=false;
    bool Va2=false;
    // check fir virtual atoms definitions (when defining com or cogs)
    if ( int(s1.find("va"))  >= 0)
      Va1=true;
    if ( int(s2.find("va"))  >= 0)
      Va2=true;
    
    // number of molecules
    vector<int> molecules;
    for(int i =0; i < at.size(); i++) {
      molecules.push_back(at.mol(i)+1);
      //cout<<at.mol(i)<<endl;
    }
    
    // Define vector specifiers 
    utils::VectorSpecifier vct1(sys,pbc);
    utils::VectorSpecifier vct2(sys,pbc);
    
    int nummol=1;

    // Let's say that the user wants to specify the molecules for
    // which the average has to be calculated.
    // Then, there will be two cases. If the flag @average is activated
    //specific molecules are from @atom definition defined, then the average runs only over those molecules.
    // else only the molecule defined by ax1 and ax2 are computed

    //fill Vector of Vectorspecifiers 
    if(args.count("average")>=0){
      nummol = molecules.size();
      for(unsigned int m = 0; m < molecules.size(); ++m) {
	ostringstream t1, t2;
	if (!Va1 && !Va2 ) {
	  t1 << s1.substr(0, s1.find("(") + 1) << molecules[m] << s1.substr(s1.find(":"), s1.size());
	  t2 << s2.substr(0, s2.find("(") + 1) << molecules[m] << s2.substr(s2.find(":"), s2.size());
	}
	else if  (Va1 && !Va2) {
	  t1 << s1.substr(0, s1.find("(") + 1) << molecules[m] 
	     << s1.substr(s1.find(":"),s1.find(",")-s1.find(":")+1)<<molecules[m]
	     <<s1.substr(s1.find(":",s1.find(":")+1), s1.size());
	  t2 << s2.substr(0, s2.find("(") + 1) << molecules[m] << s2.substr(s2.find(":"), s2.size());
	}
	else if (!Va1 && Va2) {
	  
	  t1 << s1.substr(0, s1.find("(") + 1) << molecules[m] << s1.substr(s1.find(":"), s1.size());
	  t2 << s2.substr(0, s2.find("(") + 1) << molecules[m] 
	     << s2.substr(s2.find(":"),s2.find(",")-s2.find(":")+1)<<molecules[m]
	     <<s2.substr(s2.find(":",s2.find(":")+1), s2.size());
	}
	else {
	  t1 << s1.substr(0, s1.find("(") + 1) << molecules[m] 
	     << s1.substr(s1.find(":"),s1.find(",")-s1.find(":")+1)<<molecules[m]
	     <<s1.substr(s1.find(":",s1.find(":")+1), s1.size());
	  t2 << s2.substr(0, s2.find("(") + 1) << molecules[m] 
	     << s2.substr(s2.find(":"),s2.find(",")-s2.find(":")+1)<<molecules[m]
	     <<s2.substr(s2.find(":",s2.find(":")+1), s2.size());
	}
	vct1.setSpecifier(t1.str());
	vs1.push_back(vct1);
	vct2.setSpecifier(t2.str());
	vs2.push_back(vct2);
	
      }
    } 
    else {
      vct1.setSpecifier(s1);
      vct2.setSpecifier(s2);
      vs1.push_back(vct1);
      vs2.push_back(vct2);
    }
    
    cout<<nummol<<endl;
    // prepare vector to store all data.
    double rad2= pore_radius*pore_radius;
    double  max_z,min_z,center_x,center_y;
    vector<vector<Vec> > data(nummol*3);
    // 2D vector with all selected atoms tracking if atoms were inside (1) or outside the pore (0)
    vector<vector <bool> > in_out(nummol);
    vector<vector<int > >  loads_in_time;
    
    // define input coordinate
    InG96 ic;

    int numFrames=0;
    int frames=0;
    int time_counter=0;
    int init_counter=0;

    // loop over all trajectories
    cout<<"#computing..."<<endl;
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys >> time;
	time_counter++;
	if (init_counter < first)
	  init_counter++;
	if (time_counter >= first && time_counter <= last ) {
	  frames++;
	  (*pbc.*gathmethod)();
	   // collect data before fitting
	  for(int m=0; m<nummol; m++){
	  //calculate the vectors for this molecule
	    Vec v1 = vs1[m]();
	    Vec v2 = vs2[m]();
	    //normalize these vectos
	    v1 = v1.normalize();
	    v2 = v2.normalize();
	    Vec v3 = v1.cross(v2);
	    data[3*m  ].push_back(v1);
	    data[3*m+1].push_back(v2);
	    data[3*m+2].push_back(v3);
	  }
	  // fit with respect to reference atoms
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
	  center_x /= pore_atoms.size();
	  center_y /= pore_atoms.size();
	  // loop through selected atoms to get load indexes
	  vector<int> load_indexes;
	  //loop over the molecules
	  for(int m=0; m<nummol; m++){
	    // save loaded particles at each frame
	    if (((at.pos(m)[2] >= (min_z)) && (at.pos(m)[2] <= (max_z) ))  &&  
	    	( (at.pos(m)[0]-center_x)*(at.pos(m)[0]-center_x) + 
		  (at.pos(m)[1]-center_y)*(at.pos(m)[1]-center_y) <=  rad2) ) {
	      load_indexes.push_back(m);
	      in_out[m].push_back(true);
	    } else
	      in_out[m].push_back(false);
	  }
	  // cout<<load_indexes.size()<<" "<<endl;  
	  loads_in_time.push_back(load_indexes);
	}
      }
      ic.close();
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }


    

    //Compute autocorrelation functions


    //3d vectors to store autocorrelation functions
    vector<vector<float> > P1_lists(3, vector<float>(window,0.0));
    vector<vector<float> > P2_lists(3, vector<float>(window,0.0));
    // 1d vector to count for averaging
    vector<int> average_counter(window,0);
    
    
    //multiple origins and windowed only inside tube 
 
    if (IN)
#ifdef OMP
#pragma omp parallel for  num_threads(threads)
#endif
      for(int it = 0; it < frames-window; it++) {
	int window_counter = 0;
	for(int j = it; j < it+window; j++) { //loop through frames starting from "it"
	  for (int h =0; h < loads_in_time[it].size() ; h++) {
	    int i=loads_in_time[it][h]; // select particle h loaded at time "it"
	    if (in_out[i][j]) { // if particle i  is inside at time "j"
	      // loop trough 3 vectors
	      for(int k = 0; k < 3; k++) {
		const double inp = data[3*i+k][j].dot(data[3*i+k][it]);
		P1_lists[k][window_counter]+=inp; // 1st legendre polynomial
		P2_lists[k][window_counter]+=0.5*(3.0*inp*inp-1.0); // 2nd legendre polynomial
	      }
	      average_counter[window_counter]++;
	    }
	  }
	  window_counter++;
	}
      }
    

   //multiple origins and windowed everywhere
    if (!IN)
#ifdef OMP
#pragma omp parallel for  num_threads(threads)
#endif

    for(int it = 0; it < frames-window; it++) {
      int window_counter = 0;
      for(int j = it; j < it+window; j++) { //loop through frames starting from "it"
	for (int h =0; h < nummol ; h++) {
	  if (!in_out[h][j]) { // if particle i outside pore
	    // loop trough 3 vectors
	    for(int k = 0; k < 3; k++) {
	      const double inp = data[3*h+k][j].dot(data[3*h+k][it]);
	      P1_lists[k][window_counter]+=inp; // 1st legendre polynomial
	      P2_lists[k][window_counter]+=0.5*(3.0*inp*inp-1.0); // 2nd legendre polynomial
	    }
	    average_counter[window_counter]++;
	  }
	}
	window_counter++;
      }
    }




    // print out a header
    cout << "# Calculating the autocorrelation function of the legendre"
	 << " polynomials" << endl;
    if (IN)
      cout << "#Particles inside cylindrical  pore defined by atoms " <<atoms_count<< endl;
    cout << "# of the dot-product of three molecular axes" << endl;
    cout << "# ax1 is defined as \"" << vs1[0].toString() << "\"\n";
    cout << "# ax2 is defined as \"" << vs2[0].toString() << "\"\n";
    cout << "# ax3 is the cross product of these two" << endl;
    cout << "#" << endl;
    if(nummol>1)
      cout << "# All vectors are calculated for " << vs1.size() << " molecules\n"
	   << "# autocorrelation functions are averaged over the molecules\n#\n";
    
    cout << "#" << setw(9) << "time"
	 << setw(14) << "<p_1(ax1)>"
	 << setw(14) << "<p_2(ax1)>"
	 << setw(14) << "<p_1(ax2)>"
	 << setw(14) << "<p_2(ax2)>"
	 << setw(14) << "<p_1(ax3)>"
	 << setw(14) << "<p_2(ax3)>"
	 << endl;
    



    
    for(int i = 0; i < P1_lists[0].size(); i++) {
      cout<<setw(5)<<i*time.dt();
      for(int k=0; k<3; k++){
	const double ave1 = P1_lists[k][i]/average_counter[i];
        const double ave2 = P2_lists[k][i]/average_counter[i];
	cout<< setw(14) << ave1
	     << setw(14) << ave2;
      }
      cout<<endl;
    }

/*
    // Now we define a counter for the correct normalization
    int count = 0;

    // now calculate all the autocorrelation functions.
    for(int it=0; it< numFrames; it++){
      double frame_sum1[3]={0.0,0.0,0.0};
      double frame_sum2[3]={0.0,0.0,0.0};
      count = 0;
      #ifdef OMP
      #pragma omp parallel for
      #endif
      for(int j=0; j < numFrames-it; j++){
        double sum1[3]={0.0,0.0,0.0};
        double sum2[3]={0.0,0.0,0.0};
        #ifdef OMP
        #pragma omp critical
        #endif
        {
          count++;
        }
        for(int m=0; m < nummol; m++){
	  for(int k=0; k<3; k++){
	    const double inp = data[3*m+k][j].dot(data[3*m+k][j+it]);
	    sum1[k]+=inp;
	    sum2[k]+=0.5*(3.0*inp*inp-1.0);
           
	  }
	}
        #ifdef OMP
        #pragma omp critical
        #endif
        {
          for(int k=0; k<3; k++) {
            frame_sum1[k] += sum1[k];
            frame_sum2[k] += sum2[k];
          }
        }
      }
      // now print out the information
      cout << times[it];
      for(int k=0; k<3; k++){
        const double ave1 = frame_sum1[k] / count / nummol;
        const double ave2 = frame_sum2[k] / count / nummol;
	cout << setw(14) << ave1
	     << setw(14) << ave2;
      }
      cout << endl;
    }
*/    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

