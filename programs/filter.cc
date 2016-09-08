/**
 * @file filter.cc
 * Filters a coordinate trajectory for a specific set of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor filter
 * @section filter Filters a coordinate trajectory for a specific set of atoms
 * @author @ref co
 * @date 11-6-07
 *
 * Program filter reduces coordinate trajectory file(s) and writes out a
 * trajectory (in gromos or pdb format) or position restraints file in which
 * for every frame, the coordinates are only kept for atoms that are within a 
 * specific distance of a specified part of the system. To determine if 
 * interatomic distances are within the specified cut-off, either an atomic or a
 * charge-group based cut-off scheme can be employed. Additionally, parts of the
 * system can be specified for which in all cases, the atomic coordinates 
 * should either be kept or rejected.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider as reference part of the system&gt; </td></tr>
 * <tr><td> \@cutoff</td><td>&lt;cut-off distance (nm, default: 0.0)&gt; </td></tr>
 * <tr><td> \@select</td><td>&lt;@ref AtomSpecifier "atoms" to keep&gt; </td></tr>
 * <tr><td> \@reject</td><td>&lt;@ref AtomSpecifier "atoms" to discard&gt; </td></tr>
 * <tr><td> \@pairlist</td><td>&lt;cut-off scheme (ATOMIC (default) or CHARGEGROUP)&gt; </td></tr>
 * <tr><td> [\@notimeblock</td><td>&lt;do not write timestep block&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@outformat</td><td>&lt;@ref args::OutformatParser "output coordinates format"&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input trajectory file(s)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  filter
    @topo      ex.top
    @pbc       r
    @atoms     1:1
    @cutoff    0.5
    @select    1:1-20
    @reject    1:21-51
    @pairlist  ATOMIC
    [@notimeblock]
    [@time     0 2]
    @outformat pdb
    @traj      ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace std;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "cutoff" << "atoms" << "select"
         << "reject" << "pairlist" << "outformat" << "notimeblock" << "time";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gather method>]\n";
  usage += "\t[@atoms     <atoms to consider as reference point>]\n";
  usage += "\t[@cutoff    <cutoff distance>]\n";
  usage += "\t[@select    <atoms to keep>]\n";
  usage += "\t[@reject    <atoms not to keep>]\n";
  usage += "\t[@pairlist  <ATOMIC or CHARGEGROUP>]\n";
  usage += "\t@outformat  <output coordinates format>\n";
  usage += "\t[@notimeblock <do not write timestep block>]\n";
  usage += "\t[@time      <time and dt>]\n";  
  usage += "\t@traj       <input trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    //System osys(it.system());

    System refSys(it.system());
    
    // create Time object
    utils::Time time(args);

    // get @notimeblock argument
    bool notimeblock = false;
    if (args.count("notimeblock") >= 0)
      notimeblock = true;

    // get simulation time either from the user or from the files
    bool usertime=false;
    
    if (args.count("time") > 0) {
      usertime=true;
    }

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // define input coordinates
    InG96 ic;

    // read in the atom list the has to be kept definitely
    utils::AtomSpecifier ls(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("select"),
              to = args.upper_bound("select");
      for (; iter != to; ++iter) {
        ls.addSpecifier(iter->second);
      }
    }
    // read in the atom list for atoms to throw away for sure
    utils::AtomSpecifier rej(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("reject"),
              to = args.upper_bound("reject");
      for (; iter != to; ++iter) {
        rej.addSpecifier(iter->second);
      }
    }
    // check that there are no doubles in ls and rej
    bool warn = false;
    for (int i = 0; i < rej.size(); i++) {
      if (ls.findAtom(rej.mol(i), rej.atom(i)) != -1) {
        warn = true;
      }
    }
    if (warn) {
      cerr << "select and reject show overlap. Atoms will be rejected." << endl;
    }

    // read in the reference atoms for the additional distance criterion
    utils::AtomSpecifier ref(sys);
    double cut = 0;
    {
      Arguments::const_iterator iter = args.lower_bound("atoms"),
              to = args.upper_bound("atoms");
      for (; iter != to; ++iter) {
        ref.addSpecifier(iter->second);
      }
    }
    if (ref.size()) {
      // then we need a cutoff
      if (args.count("cutoff") <= 0)
        throw gromos::Exception("filter", "If you specify reference atoms, then you need a cutoff");
      cut = args.getValue<double>("cutoff");
    }

    // read in the type
    std::string t = "ATOMIC";
    if (args.count("pairlist") > 0) {
      if (args["pairlist"] == "CHARGEGROUP") t = args["pairlist"];
      else if (args["pairlist"] != "ATOMIC") throw gromos::Exception("filter",
              "only ATOMIC and CHARGEGROUP are accepted for pairlist");
    }

    string ext;
    OutCoordinates & oc = *OutformatParser::parse(args, ext);
    oc.open(cout);

    std::ostringstream title;

    title << "Filtered trajectory. Keeping" << endl;
    if (ls.size()) {
      vector<string> s = ls.toString();
      title << "* atoms ";
      for (unsigned int i = 0; i < s.size(); i++) {
        title << s[i] << " ";
      }
      title << endl;
    }
    if (ref.size()) {
      vector<string> s = ref.toString();
      title << "* atoms within " << cut << " of atoms ";
      for (unsigned int i = 0; i < s.size(); i++) {
        title << s[i] << " ";
      }
      title << endl;
      title << "  using a";
      if (t == "ATOMIC") {
        title << "n atomic";
      } else {
        title << " chargegroup based";
      }
      title << " cutoff" << endl;
    }
    if (rej.size()) {
      vector<string> s = rej.toString();
      title << "* as long as the above atoms do not belong to ";
      for (unsigned int i = 0; i < s.size(); i++) {
        title << s[i] << " ";
      }
      title << endl;
    }

    oc.writeTitle(title.str());

    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {

      ic.open(iter->second);
      ic.select("ALL");

      // loop over all frames
      while (!ic.eof()) {
        if (!notimeblock) {
          ic >> sys >> time;
        } else {
          ic >> sys;
        }
        (*pbc.*gathmethod)();

        if (!notimeblock) {
          oc.writeTimestep(time.steps(), time.time());
        }
        utils::AtomSpecifier rls = ls;

        // add all atoms that need to be added according to the 
        // distances to reference atoms
#ifdef OMP
#pragma omp parallel for 
#endif
        for (int i = 0; i < ref.size(); i++) {
          utils::SimplePairlist spl(sys, *pbc, cut);
          spl.setAtom(*ref.atom()[i]);
          spl.setType(t);
          spl.calc();

          if ((*ref.atom()[i]).type() != utils::spec_virtual) {
            spl.addAtom(ref.mol(i), ref.atom(i));
          }
#ifdef OMP
#pragma omp critical 
#endif
          {
            rls = rls + spl;
          }
        }

        // remove atoms that are to be rejected
        for (int i = 0; i < rej.size(); i++) {
          rls.removeAtom(rej.mol(i), rej.atom(i));
        }
        rls.sort();
        oc << rls;
      }
      ic.close();
    }
  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
