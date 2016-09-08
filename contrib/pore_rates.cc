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
#include "../src/gmath/Stat.h"

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


int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "pore" << "atoms"  
         << "radius" << "traj" << "ref" << "offset" << "threads"
	 << "first"<<"last"<<"width"<<"MolAxialLength"<<"bins"<<"axialOff";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@pore    <pore atoms>\n";
  usage += "\t@atoms  <atoms to count>\n";
  usage += "\t@radius  <pore radius>\n";
  usage += "\t@ref  <reference structure, with pore aligned in z>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@offset   <radius of region around upper and lower mouth of pore>\n";
  usage += "\t@threads  <number of threads for openmp>\n";
  usage += "\t@first   <first snapshot to consider  >\n";
  usage += "\t@last   <final snapshot to consider  >\n";
  usage += "\t@width   <width of mixing shell surroundig pore mouths   >\n";
  usage += "\t@MolAxialLength <average load for single file load (effective axial length)>\n";
  usage += "\t@bins <number of bins for passage time distributions   >\n";
  usage += "\t@axialOff <pore axis offset for permeation events counting >\n";
  

  try {
    Arguments args(argc, argv, knowns, usage);
    // get the @time argument
    utils::Time time(args);
    // get the @radius argument
    double rad=args.getValue<double>("radius");
    double rad2=rad*rad;
    
    
    int threads;
    if (args.count("threads") > 0) 
      threads=args.getValue<int>("threads");
    else
      threads=1;

    // get the @offset argument
    double axialOff=0.0;
    if ( args.count("axialOff") > 0)
      axialOff=args.getValue<double>("axialOff");
    

    // get the @axialOff argument
    double offset=0.0;
    if ( args.count("offset") > 0)
      offset=args.getValue<double>("offset");
    double offset2=offset*offset;
    
    // get the @mix_width argument
    double width=0.0;
    if ( args.count("width") > 0)
      width=args.getValue<double>("width");
    double mix_radius2=(offset+width)*(offset+width);
    
     // get the @first argument
    int first = args.getValue<int>("first");
    
    // get the @last argument
    int last = args.getValue<int>("last");
    
    int bins=args.getValue<int>("bins");
    
    // get the @last argument
    double at_axial_length = args.getValue<double>("MolAxialLength");
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
      throw gromos::Exception("pore_rates",
			      "No atoms to calculate pore loadings!");
    vector<vector<int> > states(at.size());
    vector<vector<int> > states_in(at.size());

    enum states_enum {OUT=0, UP=1, DOWN=2, IN=3, MIX_UP=4, MIX_DOWN=5,ELSE=6,UU=7,DD=8};

    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

       
    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    int frames=0;
    int time_counter=0;
    int init_counter=0;
    //cout<<"computing..."<<endl;
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
	  // do the fitting (it also center sytems on ref molecule)
	  (*pbc.*gathmethod)();
	  rf.fit(&sys);        
	  //loop over pore molecules and get max and min z
	  // init min_z and max_z to the max and min values of the box
	  double min_z = sys.box().Z();
	  double max_z = -sys.box().Z();
	  double center_x=0;
	  double center_y=0;
	  double center_z=0;
	  

	// Loop trough all pore atoms to get min and max z and center
	  //#ifdef OMP
	  //#pragma omp parallel for  num_threads(threads)	\
	  //reduction(+: center_x,center_y,center_z)
	  //#endif
	  for(int i = 0; i < pore_atoms.size(); i++) {
	    center_x+=pore_atoms.pos(i)[0];
	    center_y+=pore_atoms.pos(i)[1];
	    center_z+=pore_atoms.pos(i)[2];
	    if(pore_atoms.pos(i)[2] > max_z){
	      max_z = pore_atoms.pos(i)[2];
	    }
	    if(pore_atoms.pos(i)[2] < min_z){
	      min_z = pore_atoms.pos(i)[2];
	    }
	  }
	  
	  // add axial offset
	  min_z-=axialOff;
	  max_z+=axialOff;

	  
	  center_x/=pore_atoms.size();
	  center_y/=pore_atoms.size();
	  center_z/=pore_atoms.size();
	  
	  // loop over atoms to count and check if they are within the pore
	  Vec x_y(0.0, 0.0, 0.0);
	  Vec x_y_z_up(0.0, 0.0, 0.0);
	  Vec x_y_z_down(0.0, 0.0, 0.0);
	  
	  int in=0;
	  int up=0;
	  int mix_up=0;
	  int down=0;
	  int mix_down=0;
	  int out=0;
#ifdef OMP
#pragma omp parallel for  num_threads(threads)				\
  reduction(+: in,up,down,mix_up,mix_down,out) private (x_y,x_y_z_up,x_y_z_down)
#endif
	  
	  for(int i = 0; i < at.size(); i++) {
	    x_y[0]=at.pos(i)[0]-center_x; 
	    x_y[1]=at.pos(i)[1]-center_y;
	    x_y_z_up[0]=at.pos(i)[0]-center_x; 
	    x_y_z_up[1]=at.pos(i)[1]-center_y;
	    x_y_z_up[2]=at.pos(i)[2]-max_z;
	    x_y_z_down[0]=at.pos(i)[0]-center_x; 
	    x_y_z_down[1]=at.pos(i)[1]-center_y;
	    x_y_z_down[2]=at.pos(i)[2]-min_z;
	    //Inside pore and just above and below pore center
	    if ((at.pos(i)[2] >= (min_z)) && (at.pos(i)[2] <= (max_z) ) &&  (x_y.abs2() <=  rad2) ) {
	      states[i].push_back(IN);
	      in++;
	      if (at.pos(i)[2] > center_z && at.pos(i)[2] <= (center_z + at_axial_length ) ) {
		states_in[i].push_back(UP);
		//cout<<frames<<" UP "<<states_in[i][frames-1]<<" "<<at.mol(i)+1<<endl;
	      }
	      else if (at.pos(i)[2]< center_z && at.pos(i)[2] >= (center_z - at_axial_length ) ) {
		states_in[i].push_back(DOWN);
		//cout<<frames<<" Down "<<states_in[i][frames-1]<<" "<<at.mol(i)+1<<endl;
	      }
	      else {
		states_in[i].push_back(IN);
		//cout<<frames<<" IN "<<states_in[i][frames-1]<<" "<<at.mol(i)+1<<endl;
	      }
	    } else if (  (at.pos(i)[2] > max_z) &&  (x_y_z_up.abs2() <=  offset2) ) {
	      states[i].push_back(UP);
	      states_in[i].push_back(OUT);
	      up++;
	      //Lower  mouth
	    } else if  ( (at.pos(i)[2] < min_z )  && (x_y_z_down.abs2() <= offset2) ) {
	      states[i].push_back(DOWN);
	      states_in[i].push_back(OUT);
	      down++;
	      //MIX UP
	    } else if  ( (at.pos(i)[2] > max_z )  &&  (x_y_z_up.abs2() > offset2) && 
			 (x_y_z_up.abs2() <= mix_radius2 ) ) {
	      states[i].push_back(MIX_UP);
	      states_in[i].push_back(OUT);
	      mix_up++;
	      // MIX DOWN
	    } else if  ( (at.pos(i)[2] < min_z )  &&  (x_y_z_down.abs2() > offset2) && 
			 (x_y_z_down.abs2() <= mix_radius2 ) ) {
	      states[i].push_back(MIX_DOWN);
	      states_in[i].push_back(OUT);
	      mix_down++;
	      // OUT
	    } else {
	      states[i].push_back(OUT);
	      states_in[i].push_back(OUT);
	      out++;
	    }
	  }
	  //cout<<"FRAME "<<frames + first<<" IN "<<in<<" UP "<<up<<" DOWN "<<down<<" OUT "<<out
	  //    <<" MIX_UP "<<mix_up<<" MIX_DOWN "<<mix_down<<" Total "<<up+down+in+out+mix_up+mix_down<<endl;
	}
      }
    }
    if (init_counter < first ) {
      cout<<"#not enough snapshots, check @first argument"<<endl;
      std::exit(EXIT_FAILURE); 
    }
    
    int P_u_d =0;
    int P_d_u=0;
    map <string, Stat<double> > rates_map; //this is the way to do it

    // init each member of rates_map
    Stat<double> dummy_stat;
    rates_map["T_OUT_UD"]  = dummy_stat;
    rates_map["T_UD_OUT"]  = dummy_stat;
    rates_map["T_UD_IN"]   = dummy_stat;
    rates_map["T_perm"]    = dummy_stat;
    rates_map["T_Up_Down"] = dummy_stat;
    rates_map["T_Down_Up"] = dummy_stat;

    
    // 2D (Matrix)( vector that holds cumulative permeation events in time for up and down
    // [k][0]=UP->IN->DOWN permeation events in frame k
    // [k][1]=DOWN->IN->UP permeation events in frame k
    // [k][2]= Total permeation events in frame k
    vector<vector<int> >  permeation(frames,vector<int>(3,0));

    for(int i = 0; i < at.size(); i++) {   //loop trough atoms
      int old_state = ELSE;
      int old_state_in = ELSE;
      for(int j = 0; j < frames; j++) {   // loop trough times
	//Collect T_OUT_UD (rates from out --> UP/DOWN)
	if ( states[i][j] == OUT    &&  old_state != OUT && 
	     old_state    != MIX_UP &&  old_state != MIX_DOWN  ) 
	  for(int k = j+1; k < frames; k++) {
	    if ( states[i][k] == UP || states[i][k] == DOWN ) { 
	      //T_OUT_UD_ts
	      rates_map["T_OUT_UD"].addval(k-j);
	      break;
	    }
	  }
	//Collect T_UD_IN  or T_UD_OUT  times ( from UP/DOWN --> IN or UP/DOWN --> out  )
	if ( (states[i][j] == UP   &&  old_state != UP   && old_state != MIX_UP) || 
	     (states[i][j] == DOWN &&  old_state != DOWN && old_state != MIX_DOWN)  )
	  for(int k = j+1; k < frames; k++) {
	    if (states[i][k] == OUT ) {
	      //T_UD_OUT
	      rates_map["T_UD_OUT"].addval(k-j);
	      break;
	    }
	    if (states[i][k] == IN) {
	      // T_UD_IN_ts 
	      rates_map["T_UD_IN"].addval(k-j);
	      break;
	    }
	  }
	//Collect permeation times (rates from UP --> IN --> DOWN )
	if (states[i][j] == IN  && old_state == UP ) {
	  int old_state_Perm = old_state;
	  for(int k = j+1; k < frames; k++) {
	    if (states[i][k] == UP)
	      break;
	    if (states[i][k] == DOWN && old_state_Perm == IN ) {
	      cout<<i<<"UP-IN-DOWN from "<<j<<"to "<<k<<endl;
	      permeation[k][0]++;
	      permeation[k][2]++;
	      // Tperm_ts 
	      rates_map["T_perm"].addval(k-j);
	      break;
	    }
	    old_state_Perm=states[i][k];
	  }
	}
	//Collect permeation times (rates from DOWN --> IN --> UP )
	if (states[i][j] == IN  && old_state == DOWN ) {
	  int old_state_in = old_state;
	  for(int k = j+1; k < frames; k++) {
	    if (states[i][k] == DOWN )
	      break;
	    if ( states[i][k] == UP && old_state_in == IN) {
	      //cout<<i<< "DOWN-IN-UP from "<<j<<"to "<<k<<endl;
	      permeation[k][1]++;
	      permeation[k][2]++;
	      //Tperm_ts;
	      rates_map["T_perm"].addval(k-j);
	      break;
	    }
	    old_state_in=states[i][k];
	  }
	}
	// Collect hoping rates just at the pore center UP/DOWN
	//cout << "new "<<states_in[i][j]<<"old "<<old_state_in<<endl;
	if (states_in[i][j] == UP  && old_state_in != UP )
	  for(int k = j+1; k < frames; k++) {
	    if (states_in[i][k] == OUT)
	      break;
	    if ( states_in[i][k] == DOWN) {
	      //cout<<"T_Up_Down from/at atom"<<j<<" "<<k<<" "<<at.mol(i)+1<<endl;
	      rates_map["T_Up_Down"].addval(k-j);
	      break;
	    }
	  }
	// Collect hoping rates just at the pore center DOWN/UP
	if (states_in[i][j] == DOWN  && old_state_in != DOWN )
	  for(int k = j+1; k < frames; k++) {
	    if (states_in[i][k] == OUT ) 
	      break;
	    if ( states_in[i][k] == UP) {
	      //cout<<"T_Down_Up from/at atom"<<j<<" "<<k<<" "<<at.mol(i)+1<<endl;
	      rates_map["T_Down_Up"].addval(k-j);
	      break;
	    }
	  }
	old_state =states[i][j];
	old_state_in=states_in[i][j];
      }
    }

    ofstream output("Permeation.dat");
    
    output<<"# Cumulative pemeation events for atoms "<<atoms_count<<endl;
    output<<setw(11) <<"# Time"
	  << ' ' << setw(11) <<"Up Down"
	  << ' ' << setw(11) <<"Down Up"
	  << ' ' << setw(11) <<"Total"
	  << "\n";
    output<<setw(11) <<0
	  << ' ' << setw(11) <<"0"
	  << ' ' << setw(11) <<"0"
	  << ' ' << setw(11) <<"0"
	  << "\n";
    
    output.precision(3);
    int cumulative_up_down=0;
    int cumulative_down_up=0;
    int cumulative_total=0;
    for (int i=0; i< frames; i++) {
      cumulative_up_down+=permeation[i][0];
      cumulative_down_up+=permeation[i][1];
      cumulative_total+=permeation[i][2];
      output<<setw(11) <<i+1+first
	    << ' ' << setw(11) <<cumulative_up_down
	    << ' ' << setw(11) <<cumulative_down_up
	    << ' ' << setw(11) <<cumulative_total
	    << "\n";
    }
    output.close();

    // Final rates

    cout<<"\n#Average passage times for atoms "<<atoms_count<<endl;
    cout<<setw(16) <<"# Tau"
	<< ' ' << setw(7) <<"n"
	<< ' ' << setw(11) <<"maxValAt"
	<< ' ' << setw(11) <<"average"
	<< ' ' << setw(11) <<"rmsd"
	<< ' ' << setw(11) <<"error"
	<< ' ' << setw(11) <<"min"
	<< ' ' << setw(11) <<"max"
	<< "\n";


    std::map<string,Stat<double> >::iterator data_iterator;
    for (std::map<string,Stat<double> >::iterator data_iterator=rates_map.begin(); 
	 data_iterator!=rates_map.end(); ++data_iterator) {
      if (!data_iterator->second.n()) {
	cout << setw(16) << data_iterator->first
	     << ' ' << setw(7) <<  data_iterator->second.n()
	     <<"\n";
      }
      else {
	string name = data_iterator->first +".dat"; 
	ofstream output2(name.c_str());
	output2<<"# Normalized Distribution of  "
	       << data_iterator->first<<"for atom "<<atoms_count<<endl;
	double min = data_iterator->second.min();
	double max = data_iterator->second.max();
	data_iterator->second.dist_init(data_iterator->second.min(),
					data_iterator->second.max(),bins);
	(data_iterator->second.distribution()).write_normalized(output2);
	output2.close();
	cout << setw(16) << data_iterator->first
	     << ' ' << setw(7)  << data_iterator->second.n()
	     << ' ' << setw(11) << (data_iterator->second.distribution()).maxValAt()
	     << ' ' << setw(11) << data_iterator->second.ave()
	     << ' ' << setw(11) << data_iterator->second.rmsd()
	     << ' ' << setw(11) << data_iterator->second.ee()
	     << ' ' << setw(11) << data_iterator->second.min()
	     << ' ' << setw(11) << data_iterator->second.max()
	     << "\n";
      }
    }
    cout<<"end"<<endl;
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


