// gio_OutPdb.cc

#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#include "OutPdb.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"
#include "../gcore/Box.h"
#include "../gcore/Bond.h"
#include "../gcore/Constraint.h"
#include "../utils/AtomSpecifier.h"

using gio::OutPdb;
using gio::OutPdb_i;
using namespace gcore;
using namespace std;
using namespace utils;

class OutPdb_i {
  friend class gio::OutPdb;
  ostream &d_os;
  int d_count, d_resoff, d_switch;
  double d_factor;

  OutPdb_i(ostream &os) :
  d_os(os), d_count(0), d_switch(0), d_factor(10.0) {
  }

  ~OutPdb_i() {
  }

  void writeSingleM(const Molecule &mol, const int mn);
  void writeSingleS(const Solvent &sol);
  void writeConect(const gcore::System &sys);
  void writeCryst(const gcore::Box &box);
  void writeAtomSpecifier(const AtomSpecifier & atoms);
};

OutPdb::OutPdb(ostream &os, double factor) :
OutCoordinates(), d_this(new OutPdb_i(os)), factor(factor) {
  d_this->d_factor = factor;
}

OutPdb::OutPdb(double factor) :
OutCoordinates(), factor(factor) {
  d_this = 0;
}

OutPdb::~OutPdb() {
  if (d_this)delete d_this;
}

void OutPdb::writeTitle(const string &title) {
  istringstream t(title);
  string line;
  while(!t.eof()) {
    getline(t, line);
    if (line[line.size()-1] == '\n')
      line[line.size()-1] = '\0';
    d_this->d_os << "TITLE " << line << "\n";
  }
}

void OutPdb::writeTimestep(const int step, const double time) {
  d_this->d_os.precision(9);
  d_this->d_os.setf(std::ios::fixed, std::ios::floatfield);

  d_this->d_os << "REMARK   1  TIMESTEP\t"
          << std::setw(18)
          << step
          << std::setw(20)
          << time
          << "\n";
}

void OutPdb::select(const string &thing) {
  if (thing == "ALL") {
    d_this->d_switch = 1;
  } else if (thing == "SOLVENT") {
    d_this->d_switch = 2;
  } else {
    d_this->d_switch = 0;
  }
}

void OutPdb::open(ostream &os) {
  if (d_this) {
    delete d_this;
  }
  d_this = new OutPdb_i(os);
  d_this->d_factor = factor;
}

void OutPdb::close() {
  if (d_this) delete d_this;
  d_this = 0;
}

OutPdb &OutPdb::operator<<(const gcore::System &sys) {
  d_this->d_os << "MODEL\n";
  if (sys.hasBox)
    d_this->writeCryst(sys.box());

  d_this->d_count = 0;
  d_this->d_resoff = 1;
  if (d_this->d_switch <= 1)
    for (int i = 0; i < sys.numMolecules(); ++i)
      d_this->writeSingleM(sys.mol(i), i + 1);
  if (d_this->d_switch >= 1)
    for (int i = 0; i < sys.numSolvents(); ++i)
      d_this->writeSingleS(sys.sol(i));

  // introduce another switch to turn on and off 
  // the writing of the CONECT entries?
  // --Clara
  d_this->writeConect(sys);
  d_this->d_os << "ENDMDL\n";

  return *this;
}

OutPdb &OutPdb::operator<<(const utils::AtomSpecifier & atoms) {
  d_this->d_os << "MODEL\n";
  if (atoms.sys()->hasBox)
    d_this->writeCryst(atoms.sys()->box());

  d_this->writeAtomSpecifier(atoms);
  d_this->d_os << "ENDMDL\n";

  return *this;
}

void OutPdb_i::writeSingleM(const Molecule &mol, const int mn) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  for (int i = 0; i < mol.numAtoms(); ++i) {
    ++d_count;
    int res = mol.topology().resNum(i);
    d_os << "ATOM";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(7) << d_count;
    d_os.setf(ios::left, ios::adjustfield);
    if (mol.topology().atom(i).name().length() == 4) {
      d_os << " " << setw(5) << mol.topology().atom(i).name().substr(0, 4).c_str();
    } else {
      d_os << "  " << setw(4) << mol.topology().atom(i).name().substr(0, 3).c_str();
    }
    d_os << setw(4) << mol.topology().resName(res).substr(0, 4).c_str();
    char chain = ('A' + mn - 1);
    // overflow!
    if (chain < 'A' || chain > 'Z') chain = 'Z';
    d_os << setw(1) << chain;
    d_os.setf(ios::right, ios::adjustfield);
    int resn = res + d_resoff;
    if (resn > 9999) resn = 9999;
    d_os << setw(4) << resn << "    "
            << setw(8) << mol.pos(i)[0]*d_factor
            << setw(8) << mol.pos(i)[1]*d_factor
            << setw(8) << mol.pos(i)[2]*d_factor
            << "  1.00  0.00" << endl;
  }
  d_os << "TER\n";
  d_resoff += mol.topology().numRes();
}

void OutPdb_i::writeSingleS(const Solvent &sol) {
  int na = sol.topology().numAtoms();

  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  for (int i = 0; i < sol.numPos(); ++i) {
    ++d_count;
    int res = i / na;
    // DW : current solv topo knows only the atom names for the numAtoms() atoms
    int nameid = i % na;

    d_os << "ATOM";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(7) << d_count;
    d_os.setf(ios::left, ios::adjustfield);
    //if (sol.topology().atom(i).name().length() == 4) {
    if (sol.topology().atom(nameid).name().length() == 4) {
      //d_os << " " << setw(5) << sol.topology().atom(i).name().substr(0, 4).c_str();
      d_os << " " << setw(5) << sol.topology().atom(nameid).name().substr(0, 4).c_str();
    } else {
      //d_os << "  " << setw(4) << sol.topology().atom(i).name().substr(0, 3).c_str();
      d_os << "  " << setw(4) << sol.topology().atom(nameid).name().substr(0, 3).c_str();
    }
    d_os << setw(4) << sol.topology().solvName().substr(0, 4).c_str();
    d_os.setf(ios::right, ios::adjustfield);
    //d_os << setw(5) << res+d_resoff << "    "
    d_os << setw(5) << res + 1 << "    "
            << setw(8) << sol.pos(i)[0]*d_factor
            << setw(8) << sol.pos(i)[1]*d_factor
            << setw(8) << sol.pos(i)[2]*d_factor
            << "  1.00  0.00" << endl;
  }
  d_os << "TER\n";
  d_resoff += sol.numPos() / na;
}

// --Clara

void OutPdb_i::writeConect(const gcore::System &sys) {
  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      d_os << setw(6) << "CONECT"
              << setw(5) << bit()[0] + offatom
              << setw(5) << bit()[1] + offatom
              << endl;


      ++count;
    }
    offatom += sys.mol(i).numAtoms();
  }
}

void OutPdb_i::writeCryst(const gcore::Box& box) {
  if (box.ntb() == Box::vacuum || box.ntb() == Box::truncoct)
    return;

  ++d_count;
  d_os.setf(ios::unitbuf);
  d_os.setf(ios::left, ios::adjustfield);
  d_os << setw(6) << "CRYST1";
  d_os.setf(ios::fixed | ios::right);
  d_os.precision(3);
  d_os << setw(9) << box.K().abs()*d_factor<< setw(9) << box.L().abs()*d_factor << setw(9) << box.M().abs()*d_factor;
  d_os.precision(2);
  d_os << setw(7) << box.alpha() << setw(7) << box.beta() << setw(7) << box.gamma();
  d_os.setf(ios::left, ios::adjustfield);
  d_os << setw(11) << " P 1";
  d_os.setf(ios::fixed | ios::right);
  d_os << setw(4) << 1;
  d_os << endl;
}

void OutPdb_i::writeAtomSpecifier(const AtomSpecifier& atoms) {
  const System & sys = *(atoms.sys());
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  int res = 0;
  int count = 0;
  int resoff = 0;

  for (int i = 0; i < atoms.size(); ++i) {

    int maxmol = atoms.mol(i);
    if (maxmol < 0) maxmol = sys.numMolecules();
    count = atoms.atom(i);
    resoff = 0;
    for (int j = 0; j < maxmol; j++) {
      count += sys.mol(j).numAtoms();
      resoff += sys.mol(j).topology().numRes();
    }

    if (atoms.mol(i) < 0) res = atoms.atom(i) / sys.sol(0).topology().numAtoms();
    else res = sys.mol(atoms.mol(i)).topology().resNum(atoms.atom(i));
    d_os << "ATOM";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(7) << count + 1;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << "  " << setw(4) << atoms.name(i).substr(0, 3).c_str();
    if (atoms.mol(i) < 0) d_os << setw(4) << "SOLV";
    else d_os << setw(4) << sys.mol(atoms.mol(i)).topology().resName(res).substr(0, 4).c_str();
    d_os.setf(ios::right, ios::adjustfield);

    int resn = res + resoff + 1;
    if (resn > 9999) resn = 9999;

    d_os << " " << setw(4) << resn << "    "
            << setw(8) << atoms.pos(i)[0]*d_factor
            << setw(8) << atoms.pos(i)[1]*d_factor
            << setw(8) << atoms.pos(i)[2]*d_factor
            << "  1.00  0.00" << endl;
  }
  // now get the bonds
  int molcount = 0;

  for (int m = 0; m < sys.numMolecules(); m++) {
    BondIterator bi(sys.mol(m).topology());
    for (; bi; ++bi) {
      int index_i = atoms.findAtom(m, bi()[0]);
      int index_j = atoms.findAtom(m, bi()[1]);
      int count_i = 0, count_j = 0;

      if (index_i != -1 && index_j != -1) {
        count_i = atoms.atom(index_i) + molcount + 1;
        count_j = atoms.atom(index_j) + molcount + 1;

        d_os << "CONECT " << count_i << " " << count_j << endl;

      }
    }
    molcount += sys.mol(m).numAtoms();

  }
  //also for solvent
  int oldoffset = -1;

  for (int j = 0; j < atoms.size(); j++) {
    // check if it is a solvent
    if (atoms.mol(j) < 0) {
      // get the atom number in this solvent
      int si = atoms.atom(j) % sys.sol(0).topology().numAtoms();
      int offset = atoms.atom(j) - si;
      if (offset != oldoffset) {
        //new molecule
        ConstraintIterator ci(sys.sol(0).topology());
        for (; ci; ++ci) {
          int index_i = atoms.findAtom(atoms.mol(j), offset + ci()[0]);
          int index_j = atoms.findAtom(atoms.mol(j), offset + ci()[1]);
          int count_i = 0, count_j = 0;

          if (index_i != -1 && index_j != -1) {
            count_i = atoms.atom(index_i) + molcount + 1;
            count_j = atoms.atom(index_j) + molcount + 1;

            d_os << "CONECT " << count_i << " " << count_j << endl;
          }
        }
      }
      oldoffset = offset;
    }
  }
  d_os << "TER\n";
}
