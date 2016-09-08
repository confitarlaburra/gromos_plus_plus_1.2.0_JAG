/**
 * @file dg_ener.cc
 * Applies the perturbation formula based on two lists of energies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dg_ener
 * @section dg_ener Applies the perturbation formula based on two lists of energies
 * @author @ref co, sr
 * @date 31-10-08 - changed 25-08-11
 *
 * Program dg_ener applies the perturbation formula to calculate the free 
 * energy difference between two states A and B. It reads in the output of 
 * program @ref ener, which can be calculated for the same trajectory using 
 * two different Hamiltonians. The free energy difference is calculated as
 * @f[ \Delta G_{BA} = -k_B T \ln < e^{-(H_B - H_A)/k_B T} > @f]
 * where the average is over all entries of the energy files that are specified
 * and the Hamiltonians are taken from the last column of these files.
 * 
 * For every line in the energy files, the program writes out the energy
 * difference and the Boltzmann probability for that particular frame. For
 * numerical reasons, this probability is normalized such that the highest
 * probability has a value of 1.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@temp</td><td>&lt;temperature for perturbation&gt; </td></tr>
 * <tr><td> \@stateA</td><td>&lt;energy file for state A&gt; </td></tr>
 * <tr><td> \@stateB</td><td>&lt;energy file for state B&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  dg_ener
    @temp    300
    @stateA  ener_output_A.dat
    @stateB  ener_output_B.dat
 @endverbatim
 *
 * <hr>
 */

#include "../src/args/Arguments.h"
#include "../src/gmath/Physics.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <sstream>

using namespace args;
using namespace std;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "temp" << "stateA" << "stateB";

  string usage = argv[0];
  usage += "\n\t@temp <temperature for perturbation>\n";
  usage += "\t@stateA <energy file for state A>\n";
  usage += "\t@stateB <energy file for state B>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // set some values
    double temp = args.getValue<double>("temp");
    ifstream stateA, stateB;

    //open files
    {
      Arguments::const_iterator iter=args.lower_bound("stateA");
      if(iter!=args.upper_bound("stateA")){
        stateA.open((iter->second).c_str());
        if(!stateA)
          throw gromos::Exception("dg_ener", "could not open energy file for state A\n");
      }
      else
        throw gromos::Exception("dg_ener", "energy file for state A missing\n");
      iter=args.lower_bound("stateB");
      if(iter!=args.upper_bound("stateB")){
	stateB.open((iter->second).c_str());
        if(!stateB)
          throw gromos::Exception("dg_ener", "could not open energy file for state B\n");
      }
      else
	throw gromos::Exception("dg_ener", "energy file for state B missing\n");
    }
    cout.precision(12);
    
    //print title
    cout << "# Time"
	 << setw(12) << "DE_tot"
	 << setw(12) << "probability"
	 << setw(12) << "DG_BA" 
	 << endl;
    
    string sdum;
    double timeA, timeB;
    
    bool errorA = false, errorB = false;
    bool eofA = false, eofB = false;
    bool timeWarning = false;

    vector<double> delta_v;
    vector<double> t;

    double covA, nbA, totA, covB, nbB, totB;

    while(true){

      // read from state A
      while(true){
	std::getline(stateA, sdum);
	if (stateA.eof()){
	  eofA = true;
	  break;
	}
	std::string::size_type it = sdum.find('#');
	if (it != std::string::npos)
	  sdum = sdum.substr(0, it);
	if (sdum.find_first_not_of(" \t") == std::string::npos)
	  continue;
	std::istringstream is(sdum);
	if (!(is >> timeA >> covA >> nbA >> totA))
	  errorA = true;
	break;
      }
      // read from state B
      while(true){
	std::getline(stateB, sdum);
	if (stateB.eof()){
	  eofB = true;
	  break;
	}
	std::string::size_type it = sdum.find('#');
	if (it != std::string::npos)
	  sdum = sdum.substr(0, it);
	if (sdum.find_first_not_of(" \t") == std::string::npos)
	  continue;
	std::istringstream is(sdum);
	if (!(is >> timeB >> covB >> nbB >> totB))
	  errorB = true;
	break;
      }
      // check eof / error
      if (eofA && eofB) break;
      if (eofA || eofB || errorA || errorB)
	throw gromos::Exception("dg_ener", "Error while reading file for state A or state B: check if number of lines or columns are identical");
      if (timeA != timeB)
	timeWarning = true;
      
      // now we have two states read
      const double dh = totB - totA;
      delta_v.push_back(dh);
      
      t.push_back(timeA);
    }

    double p, dg=0.0;
    const double beta = 1 / (gmath::physConst.get_boltzmann() * temp);
    double sum = -delta_v[0]*beta;
    
    // now loop over all values that were read in 
    for(unsigned int i=1; i<delta_v.size(); i++){

      p = -delta_v[i]*beta;
      sum = std::max(sum,p) + log(1 + exp(std::min(sum,p) - std::max(sum,p)));
      dg = - (sum - log(i)) / beta;
      
      cout.precision(5);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << t[i]
           << setw(12) << delta_v[i]
           << setw(12) << exp(p)
           << setw(12) << dg
           << endl;
    }

    cout << "# final result dG_BA = G_B - G_A: " << dg << endl;

    stateA.close();
    stateB.close();

    if (timeWarning){
      cout << "#\n# WARNING: time was not equal in state A and state B!" << endl;
      cerr << "\nWARNING: time was not equal in state A and state B!" << endl;
    }
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    return 1;
  }
  return 0;
}
