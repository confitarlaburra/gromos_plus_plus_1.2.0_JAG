#include <cassert>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;


set<int> neighbours(int a, vector<Bond> & bond);

int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "topo" << "pos" << "pbc";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@pos <coordinates>\n";
  usage += "\t@pbc <periodic boundary conditions>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());

    // read coordinates
    InG96 ic(args["pos"]);
    ic >> sys;

    // parse boundary conditions
    bound::Boundary *pbc = BoundaryParser::boundary(sys, args);
    
    // create a linear topology
    System syo;
    gcore::LinearTopology lt(syo);

    utils::AtomSpecifier as(sys);
    for(int i=0; i< sys.numMolecules(); i++){
      for(int j=0; j< sys.mol(i).numAtoms(); j++){
	as.addAtom(i,j);
      }
    }

    //std::cout << "atoms " << as.size() << std::endl;
    
    // make a double loop over all atoms to find the bonds
    vector<Bond> bonds;
    for(int i=0; i< as.size(); i++){
      for(int j=i+1; j < as.size(); j++){
	gmath::Vec v = as.pos(i) - pbc->nearestImage(as.pos(i), as.pos(j), sys.box());
	if(v.abs2() < 0.0400){
	  Bond b(i,j);
	  b.setType(14);
	  
	  bonds.push_back(b);
	  lt.addBond(b);
	  
	}
	
      }
    }

    //std::cout << "bonds " << bonds.size() << std::endl;
    //std::cout << "bonds " << lt.bonds().size() << std::endl;

    //get the exclusions (everything is aromatic) and set the atoms
    for(int i=0; i< as.size(); i++){
      AtomTopology a=sys.mol(as.mol(i)).topology().atom(as.atom(i));
      
      set<int> exclusions;
      set<int> nb1=neighbours(i, bonds);
      set<int> nb2;
      set<int> nb3;
      set<int> tmp;
      
      for(set<int>::iterator it=nb1.begin(), to=nb1.end(); it!=to; ++it){
	exclusions.insert(*it);
	nb2=neighbours(*it, bonds);
	for(set<int>::iterator it2=nb2.begin(), to2=nb2.end(); it2!=to2; ++it2){
	  exclusions.insert(*it2);
	//removed 1-4 exclusions JAG 9.11.12 (just for 5.5 tubes)
	  tmp=neighbours(*it2, bonds);
	  for(set<int>::iterator it3=tmp.begin(), to3=tmp.end();
	      it3!=to3; ++it3){
	       exclusions.insert(*it3);
	  }
	}
      }

      Exclusion ex, ex14;
      for(set<int>::iterator it=exclusions.begin(), to=exclusions.end();
	  it!=to; ++it)
	if(*it > int(i)) ex.insert(*it);
	
      a.setExclusion(ex);
      a.setExclusion14(ex14);
      
      lt.addAtom(a);
      lt.setResNum(i,0);
    }


    //cout << "did the atoms" << std::endl;

    // now the angles, and the impropers on the atoms
    for(int i=0; i< as.size(); i++){
      set<int> nb=neighbours(i,bonds);
      vector<int> nb2;
      
      for(set<int>::iterator itt=nb.begin(), to=nb.end(); itt!=to; ++itt)
	nb2.push_back(*itt);
   
      //cout << "atom " << i+1 << " has " << nb.size() << " neighbours" << std::endl;	
      if(nb.size()==3){
        Angle a1(nb2[0], i, nb2[1]);
        Angle a2(nb2[0], i, nb2[2]);
        Angle a3(nb2[1], i, nb2[2]);
      
        a1.setType(25);
        a2.setType(25);
        a3.setType(25);
        lt.addAngle(a1);
        lt.addAngle(a2);
        lt.addAngle(a3);
	// if we use charmm parameters we shoul add torsions and not the impropers
        Improper ii(i, nb2[0], nb2[1], nb2[2]);
        ii.setType(0);
        lt.addImproper(ii);
      }
      else if(nb.size()==2){
        Angle a1(nb2[0], i, nb2[1]);
      
        a1.setType(25);
        lt.addAngle(a1);   
      }
    }
    // and the impropers along the bonds (have to be defined as 0 not 180!)
    for(int i=0; i<bonds.size(); i++){
      int a,a2, b,c,d, d2;
      bool skip=false;
      b = bonds[i][0];
      c = bonds[i][1];
      
      set<int> nb1=neighbours(b, bonds);
      set<int> nb2=neighbours(c, bonds);
      a=c;
      a2=c;
      d=b;
      d2=b;
      
      // select an initial guess
      set<int>::iterator it1=nb1.begin();
      set<int>::iterator it2=nb2.begin();

      while(a==c){
	a=*it1;
	++it1;
      }
      if(nb1.size()>2){
        it1=nb1.begin();
        a2=*it1;
        if(a2==a || a2==c) ++it1;
        a2=*it1;
        if(a2==a || a2==c) ++it1;
        a2=*it1;
        if(a2==a || a2==c) std::cerr << "something is wrong!" << std::endl;
        }
      else
	a2=-1;
 
      while(d==b){
	d=*it2;
	++it2;
      }
      if(nb2.size()>2){
        it2=nb2.begin();
        d2=*it2;
        if(d2==d || d2==b) ++it2;
        d2=*it2;
        if(d2==d || d2==b) ++it2;
        d2=*it2;
        if(d2==d || d2==b) std::cerr << "something is wrong!" << std::endl;
        }
      else
	d2=-1;
      // we have to calculate the current value
      
      gmath::Vec tmpA = as.pos(a) - pbc->nearestImage(as.pos(a), as.pos(b), sys.box());
      gmath::Vec tmpB = as.pos(d) - pbc->nearestImage(as.pos(d), as.pos(c), sys.box());
      gmath::Vec tmpC = as.pos(c) - pbc->nearestImage(as.pos(c), as.pos(b), sys.box());

      gmath::Vec p1 = tmpA.cross(tmpC);
      gmath::Vec p2 = tmpB.cross(tmpC);
      
      double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

      double value = acos(cosphi)*180 / M_PI;

      gmath::Vec p3 = p1.cross(p2);
      if (p2.dot(tmpA) < 0)
	value = -value;
      
      //std::cout << a << " " << b << " " << c << " " << d <<  " : " << value << std::endl;
      if(value > 90.0 || value < -90.0 ){
        if(d2!=-1)
	  d=d2;
        else if(a2!=-1)
          a=a2;
      }
      //std::cout << a << " " << b << " " << c << " " << d << std::endl;
      if(d!=-1 && a!=-1){
        Improper ii(a,b,c,d);
        ii.setType(0);
        lt.addImproper(ii);
      }
    }
    
   
    lt.setResName(0,"CCC");
    
    //std::cout << "exclusions " << std::endl;
    
    // and parse the linearized thing back into a topology (which might
    // have a zillion molecules, because bonds have been removed)
    syo = lt.parse();
    
    //std::cout << "parsed " << std::endl;
    
    // take the old solvent
    syo.addSolvent(sys.sol(0));
    
    // and write out the new topology
    OutTopology ot(cout);
    ostringstream os;
    os << "Topology based on " << args["topo"] << endl;
    os << "Added bonded interactions" << endl;
    
    ot.setTitle(os.str());
    
    ot.write(syo, it.forceField());
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}



set<int> neighbours(int a, vector<Bond> & bond)
{
  set<int> tmp;
  for(unsigned int i=0; i< bond.size(); i++){
    if(bond[i][0]==a) tmp.insert(bond[i][1]);
    if(bond[i][1]==a) tmp.insert(bond[i][0]);
  }
  return tmp;
}
