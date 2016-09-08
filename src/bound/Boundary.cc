// bound_Boundary.cc

#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <set>

#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include "../fit/PositionUtils.h"
#include "../gio/InG96.h"
#include "Boundary.h"


using gmath::Vec;
using namespace gcore;
using bound::Boundary;
using bound::Boundary_i;
using namespace std;

class Boundary_i {
  friend class bound::Boundary;
  gcore::System *d_sys;
  gcore::System *d_refSys;
  vector<const gmath::Vec *> d_ref;
  char d_type;

  Boundary_i() : d_ref() {
  }

  ~Boundary_i() {
  }
};

Boundary::Boundary(System *sys) : d_this(new Boundary_i()) {
  d_this->d_sys = sys;
  d_this->d_refSys = NULL;
  for (int i = 0; i < d_this->d_sys->numMolecules(); ++i) {
    //the system might not have coordinates yet!
    if (d_this->d_sys->mol(i).numPos() == 0)
      d_this->d_sys->mol(i).initPos();
    d_this->d_ref.push_back(&d_this->d_sys->mol(i).pos(0));
  }
}

void Boundary::setReference(int i, const Vec &vec) {
  assert(i<int(d_this->d_ref.size()));
  d_this->d_ref[i] = &vec;
}

void Boundary::setReference(System const & sys) {
  assert(int(d_this->d_ref.size()) == sys.numMolecules());
  cerr << "d_this->d_ref.size() = " << d_this->d_ref.size() << endl;
  for (int i = 0; i < sys.numMolecules(); ++i) {

    d_this->d_ref[i] = &sys.mol(i).pos(0);
  }
}

void Boundary::setReferenceFrame(std::string file) {
  gio::InG96 in(file);
  in.select("ALL");
  d_this->d_refSys = new System(sys());
  in >> refSys();
}

void Boundary::setReferenceSystem(System system) {
  if (d_this->d_refSys == NULL) {
    d_this->d_refSys = new System(system);
  } else {
    *(d_this->d_refSys) = system;
  }
}

bound::Boundary::~Boundary() {
  if (d_this->d_refSys != NULL)
    delete d_this->d_refSys;
  delete d_this;
}

const gmath::Vec &Boundary::reference(int i)const {
  assert(d_this != NULL);

  assert(i<int(d_this->d_ref.size()));
  return *d_this->d_ref[i];
}

char Boundary::type() {
  return d_this->d_type;
}

void Boundary::setType(char t) {
  switch (t) {
    case 'v':
      d_this->d_sys->box().setNtb(gcore::Box::vacuum);
      break;
    case 'r':
      d_this->d_sys->box().setNtb(gcore::Box::rectangular);
      break;
    case 't':
      d_this->d_sys->box().setNtb(gcore::Box::truncoct);
      break;
    case 'c':
      d_this->d_sys->box().setNtb(gcore::Box::triclinic);
      break;
    default:
      stringstream msg;
      msg << "periodic boundary condition '" << t << "' unknow. Known "
              "boundaries are r (rectangular), t (truncated octahedron), c (triclinic) and v (vacuum)";
      throw gromos::Exception("Boundary", msg.str());
      break;
  }
  d_this->d_type = t;

}

System &Boundary::sys() {
  return *d_this->d_sys;
}

System &Boundary::refSys() {
  return *d_this->d_refSys;
}

void Boundary::nogather() {

}

// gather based on a general list

void Boundary::gatherlist() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  // gather the first molecule
  Molecule &mol = sys().mol(0);
  mol.pos(0) = nearestImage(reference(0), mol.pos(0), sys().box());
  for (int j = 1; j < mol.numPos(); ++j)
    mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());

  // gather the rest molecules
  // check whether the molecule should be gathered according to an atom list

  for (int i = 1; i < sys().numMolecules(); ++i) {

    Molecule &mol = sys().mol(i);
    int m = int(sys().primlist[i][0]);
    int n = int(sys().primlist[i][1]);
    int o = int(sys().primlist[i][2]);

    mol.pos(m) = nearestImage(sys().mol(n).pos(o), mol.pos(m), sys().box());

    if (m == 0) {
      for (int j = 1; j < mol.numPos(); ++j)
        mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
    } else {
      for (int j = m - 1; j >= 0; --j) {
        mol.pos(j) = nearestImage(mol.pos(j + 1), mol.pos(j), sys().box());
      }
      for (int j = m + 1; j < mol.numPos(); ++j) {
        mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
      }
    }
  }

  // now calculate cog
  Vec cog(0.0, 0.0, 0.0);
  int natom = 0;
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    for (int a = 0; a < mol.numAtoms(); ++a) {
      cog += mol.pos(a);
      natom += 1;
    }
  }
  cog /= double(natom);

  // do the solvent
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(cog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
};

// gather in term of time

void Boundary::gathertime() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  if (refSys().sol(0).numPos() != sys().sol(0).numPos())
    throw gromos::Exception("Gather problem",
          "Number of solvent atoms are not equal in reference and the current system! Abort!");

  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    Molecule &refmol = refSys().mol(i);
    for (int j = 0; j < mol.numPos(); ++j) {
      mol.pos(j) = nearestImage(refmol.pos(j), mol.pos(j), sys().box());
      refmol.pos(j) = mol.pos(j);
    }
  }

  // calculate cog
  Vec cog(0., 0., 0.);
  int count = 0;
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    for (int j = 0; j < mol.numPos(); ++j) {
      cog += mol.pos(j);
      count += 1;
    }
  }
  cog /= double(count);

  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refsol = refSys().sol(0);

  if (refsol.numPos() == sol.numPos())
    for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
      sol.pos(i) = nearestImage(refsol.pos(i), sol.pos(i), sys().box());
      refsol.pos(i) = sol.pos(i);
      for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
        sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
        refsol.pos(j) = sol.pos(j);
      }
    } else {
    std::cout << "# solv num " << sol.numPos()
            << " and refsolv num " << refsol.numPos()
            << " are not equal. solv gathering based on cog : " << std::endl;

    for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
      sol.pos(i) = nearestImage(cog, sol.pos(i), sys().box());
      for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
        sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
      }
    }
  }

};

// gather the 1st frame based on an atom list, then the rest in term of time
// everytime we update the reference system, thus to avoid changing the code for an old state

void Boundary::gatherltime() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  if (sys().primlist[0][0] == 31415926) {
    for (int i = 0; i < sys().numMolecules(); ++i) {
      Molecule &mol = sys().mol(i);
      Molecule &refmol = refSys().mol(i);
      mol.pos(0) = nearestImage(refmol.pos(0), mol.pos(0), sys().box());
      refmol.pos(0) = mol.pos(0);
      for (int j = 1; j < mol.numPos(); ++j) {
        mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
        refmol.pos(j) = mol.pos(j);
      }
    }

    // calculate the cog
    Vec cog(0., 0., 0.);
    int count = 0;
    for (int i = 0; i < sys().numMolecules(); ++i) {
      Molecule &mol = sys().mol(i);
      for (int j = 0; j < mol.numPos(); ++j) {
        cog += mol.pos(j);
        count += 1;
      }
    }
    cog /= double(count);

    // do the solvent
    Solvent &sol = sys().sol(0);
    Solvent &refsol = refSys().sol(0);

    if (refsol.numPos() == sol.numPos())
      for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
        sol.pos(i) = nearestImage(refsol.pos(i), sol.pos(i), sys().box());
        refsol.pos(i) = sol.pos(i);
        for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
          sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
          refsol.pos(j) = sol.pos(j);
        }
      } else {
      std::cout << "# solv num " << sol.numPos()
              << " and refsolv num " << refsol.numPos()
              << " are not equal. solv gathering based on cog : " << endl;

      for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
        sol.pos(i) = nearestImage(cog, sol.pos(i), sys().box());
        for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
          sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
        }
      }
    }
  } else {
    // gather the first molecule
    Molecule &mol = sys().mol(0);
    mol.pos(0) = nearestImage(reference(0), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numPos(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());

    // gather the rest molecules
    // check whether the molecule should be gathered according to an atom list

    for (int i = 1; i < sys().numMolecules(); ++i) {

      Molecule &mol = sys().mol(i);
      int m = int(sys().primlist[i][0]);
      int n = int(sys().primlist[i][1]);
      int o = int(sys().primlist[i][2]);

      mol.pos(m) = nearestImage(sys().mol(n).pos(o), mol.pos(m), sys().box());

      if (m == 0) {
        for (int j = 1; j < mol.numPos(); ++j)
          mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
      } else {
        for (int j = m - 1; j >= 0; --j) {
          mol.pos(j) = nearestImage(mol.pos(j + 1), mol.pos(j), sys().box());
        }
        for (int j = m + 1; j < mol.numPos(); ++j) {
          mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
        }
      }
    }

    // now calculate cog
    Vec cog(0.0, 0.0, 0.0);
    int natom = 0;
    for (int i = 0; i < sys().numMolecules(); ++i) {
      Molecule &mol = sys().mol(i);
      for (int a = 0; a < mol.numAtoms(); ++a) {
        cog += mol.pos(a);
        natom += 1;
      }
    }
    cog /= double(natom);

    // do the solvent
    Solvent &sol = sys().sol(0);
    for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
      //sol.pos(i)=nearestImage(reference(0),sol.pos(i),sys().box());
      sol.pos(i) = nearestImage(cog, sol.pos(i), sys().box());
      for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
        sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
      }
    }
    // Here we define the gathering of next frame wrt time rather than a list
    sys().primlist[0][0] = 31415926;
    // set the current system as reference System to gather the rest with respect to time
    setReferenceSystem(sys());
  }
};

// gather based on a reference structure

void Boundary::gatherref() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  if (sys().numMolecules() != refSys().numMolecules())
    throw gromos::Exception("Gather problem",
          "Number of SOLUTE  molecules in reference and frame are not the same.");
  if (sys().sol(0).numPos() != refSys().sol(0).numPos())
    std::cout << "# WARNING: " << endl
          << "\n# Number of SOLVENT molecules in reference and frame are not the same."
          << "\n# Gathering of solvent will be bashed on COG of solute" << std::endl;


  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    Molecule &refMol = refSys().mol(i);
    for (int j = 0; j < mol.numPos(); ++j) {
      // gather to the reference
      mol.pos(j) = nearestImage(refMol.pos(j), mol.pos(j), sys().box());
    }
  }
  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refSol = refSys().sol(0);

  Vec solcog(0., 0., 0.);
  if (sys().sol(0).numPos() == refSys().sol(0).numPos()) {
    for (int i = 0; i < refSol.numPos(); i += sol.topology().numAtoms()) {
      //sol.pos(i) = nearestImage(refSol.pos(i), sol.pos(i), sys().box());
      for (int j = i; j < (i + sol.topology().numAtoms()); ++j) {
        solcog += refSol.pos(j);
      }
    }
    solcog /= refSol.numPos();
  } else {
    int num = 0;
    for (int i = 0; i < sys().numMolecules(); ++i) {
      Molecule &mol = sys().mol(i);
      for (int j = 0; j < mol.numPos(); ++j) {
        solcog += mol.pos(j);
      }
      num += mol.numPos();
    }
    solcog /= num;
  }

  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(solcog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
};

void Boundary::gatherrtime() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  const gcore::Box & box = sys().box();

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");


  if (sys().numMolecules() != refSys().numMolecules())
    throw gromos::Exception("Gather problem",
          "Number of molecules in reference and frame are not the same.");
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    Molecule &refMol = refSys().mol(i);
    for (int j = 0; j < mol.numPos(); ++j) {
      // gather to the reference
      mol.pos(j) = nearestImage(refMol.pos(j), mol.pos(j), box);
      // set the current frame as the reference for the next frame
      refMol.pos(j) = mol.pos(j);
    }
  }
  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refSol = refSys().sol(0);
  Vec solcog(0., 0., 0.);
  if (sol.numPos() != refSol.numPos()) {
    cout << "WARNING: " << endl
            << "Number of SOLVENT molecules in reference and frame are not the same." << endl
            << "Gathering of solvent will be bashed on COG of solute" << endl;

    int num = 0;
    for (int i = 0; i < sys().numMolecules(); ++i) {
      Molecule &mol = sys().mol(i);
      for (int j = 0; j < mol.numPos(); ++j) {
        solcog += mol.pos(j);
      }
      num += mol.numPos();
    }
    solcog /= num;

    for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
      sol.pos(i) = nearestImage(solcog, sol.pos(i), sys().box());
      for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
        sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
      }
    }
  } else {
    for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
      sol.pos(i) = nearestImage(refSol.pos(i), sol.pos(i), sys().box());
      refSol.pos(i) = sol.pos(i);
      for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
        sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
        refSol.pos(j) = sol.pos(j);
      }
    }
  }
}

// gather based on bond connection

void Boundary::gatherbond() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(i), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numPos(); ++j) {

      //find a previous atom to which we are bound
      BondIterator bi(mol.topology());
      int k = 0;
      for (; bi; ++bi)
        if (bi()[1] == j) {
          k = bi()[0];
          break;
        }
      mol.pos(j) = nearestImage(mol.pos(k), mol.pos(j), sys().box());
    }
  }

  int num = 0;
  Vec cog(0., 0., 0.);
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    for (int j = 0; j < mol.numPos(); ++j) {
      cog += mol.pos(j);
    }
    num += mol.numPos();
  }
  cog /= num;

  // do the solvent
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(cog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
};

void Boundary::gather() {

  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(0), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numPos(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
  }
  // do the solvent 
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(reference(0), sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
}

void Boundary::gathergr() {

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(i), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numAtoms(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
  }
}

void Boundary::gathermgr() {

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  const Vec centre(0.5 * sys().box().K().abs(), 0.5 * sys().box().L().abs(), 0.5 * sys().box().M().abs());

  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(i), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numAtoms(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());

    // now calculate cog
    Vec cog(0.0, 0.0, 0.0);
    for (int a = 0; a < mol.numAtoms(); ++a)
      cog += mol.pos(a);
    cog /= mol.numAtoms();
    Vec cog_box = nearestImage(centre, cog, sys().box());
    Vec trans = cog_box - cog;
    fit::PositionUtils::translate(mol, trans);
  }
}

void Boundary::coggather() {

  //std::cout << "# gather with cog " << std::endl;
  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  Molecule &mol = sys().mol(0);
  Solvent &sol = sys().sol(0);

  Vec ref(0.0, 0.0, 0.0);
  Vec cog(0.0, 0.0, 0.0);
  int atoms = 0;

  // do mol(0) with respect to ref (0,0,0)
  mol.pos(0) = nearestImage(ref, mol.pos(0), sys().box());
  for (int j = 1; j < mol.numAtoms(); ++j) {
    mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
  }

  // calculate COG of mol(0)
  for (int i = 0; i < mol.numAtoms(); i++) {
    cog = cog + mol.pos(i);
    ++atoms;
  }
  cog = (1.0 / double(atoms)) * cog;

  // do the rest of the molecules
  for (int i = 1; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(cog, mol.pos(0), sys().box());
    for (int j = 1; j < mol.numPos(); ++j) {
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
    }
  }

  // do the solvent 
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(cog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }

}

void Boundary::crsgather() {

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(i), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numAtoms(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Vec cog(0.0, 0.0, 0.0);
    int numat = 0;
    for (int j = 0; j < sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0 / double(numat)) * cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nearestImage((0.0,0.0,0.0),vcog[0],sys().box());

  // Now gather cog's w.r.t. cog of previous molecule
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog = vcog[0];
  for (int i = 0; i < sys().numMolecules() - 1; ++i) {
    Vec nimj = nearestImage(vcog[i], vcog[i + 1], sys().box());
    vcog[i + 1] = nimj;
    ocog += vcog[i + 1];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    for (int j = 0; j < mol.numAtoms(); ++j) {
      mol.pos(j) = nearestImage(vcog[i], mol.pos(j), sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog = ocog / double(sys().numMolecules());
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(ocog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
}

void Boundary::seqgather() {

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(i), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numAtoms(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Vec cog(0.0, 0.0, 0.0);
    int numat = 0;
    for (int j = 0; j < sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0 / double(numat)) * cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nearestImage((0.0,0.0,0.0),vcog[0],sys().box());

  // Now gather cog's w.r.t. cog of previous molecule
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog = vcog[0];
  for (int i = 0; i < sys().numMolecules() - 1; ++i) {
    // crs:
    // Vec nimj=nearestImage(vcog[i],vcog[i+1],sys().box());
    // seq:
    Vec nimj = nearestImage(ocog / double(i + 1), vcog[i + 1], sys().box());
    vcog[i + 1] = nimj;
    ocog += vcog[i + 1];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    for (int j = 0; j < mol.numAtoms(); ++j) {
      mol.pos(j) = nearestImage(vcog[i], mol.pos(j), sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog = ocog / double(sys().numMolecules());
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(ocog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
}

void Boundary::gengather() {

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(i), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numAtoms(); ++j)
      mol.pos(j) = nearestImage(mol.pos(j - 1), mol.pos(j), sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Vec cog(0.0, 0.0, 0.0);
    int numat = 0;
    for (int j = 0; j < sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0 / double(numat)) * cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nearestImage((0.0,0.0,0.0),vcog[0],sys().box());

  // Use vcogi to make the graph connecting the closest cog's
  std::vector<int> vcogi;
  for (int i = 0; i < sys().numMolecules(); ++i) {
    vcogi.push_back(i);
  }

  // Now gather cog's w.r.t. each other
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog = vcog[0];
  for (int i = 0; i < sys().numMolecules() - 1; ++i) {
    // Determine closest nim to i among remaining molecules (using vcogi)
    int bufi = vcogi[i];
    int inimcogi = vcogi[i + 1];
    Vec nimcogi = nearestImage(vcog[bufi], vcog[inimcogi], sys().box());
    int jclose = i + 1;
    for (int j = i + 2; j < sys().numMolecules(); ++j) {
      int bufj = vcogi[j];
      if ((nearestImage(vcog[bufi], vcog[bufj], sys().box()) - vcog[bufi]).abs()<(nimcogi - vcog[bufi]).abs()) {
        nimcogi = nearestImage(vcog[bufi], vcog[bufj], sys().box());
        inimcogi = bufj;
        jclose = j;
      }
    }
    // Now swap inimcogi with i+1 in vcogi
    int bufci = vcogi[i + 1];
    vcogi[i + 1] = inimcogi;
    vcogi[jclose] = bufci;

    // Set vcog[i+1] either to its nim to vcog[i], or to
    // nim to overall cog of molecules[1 ... i], depending
    // on what corresponds with the closest distance
    Vec nic1 = nimcogi;
    Vec nic2 = nearestImage(ocog / double(i + 1), nimcogi, sys().box());
    if ((nic1 - vcog[bufi]).abs()<(nic2 - ocog / double(i + 1)).abs()) {
      vcog[inimcogi] = nic1;
    } else {
      vcog[inimcogi] = nic2;
    }
    ocog += vcog[inimcogi];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  // as a reference
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    for (int j = 0; j < mol.numAtoms(); ++j) {
      mol.pos(j) = nearestImage(vcog[i], mol.pos(j), sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog = ocog / double(sys().numMolecules());
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(ocog, sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
}

void Boundary::bondgather() {

  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");

  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    mol.pos(0) = nearestImage(reference(0), mol.pos(0), sys().box());
    for (int j = 1; j < mol.numPos(); ++j) {

      //find a previous atom to which we are bound
      BondIterator bi(mol.topology());
      int k = 0;
      for (; bi; ++bi)
        if (bi()[1] == j) {
          k = bi()[0];
          break;
        }
      mol.pos(j) = nearestImage(mol.pos(k), mol.pos(j), sys().box());
    }
  }

  // do the solvent 
  Solvent &sol = sys().sol(0);
  for (int i = 0; i < sol.numPos(); i += sol.topology().numAtoms()) {
    sol.pos(i) = nearestImage(reference(0), sol.pos(i), sys().box());
    for (int j = i + 1; j < (i + sol.topology().numAtoms()); ++j) {
      sol.pos(j) = nearestImage(sol.pos(j - 1), sol.pos(j), sys().box());
    }
  }
}

void Boundary::refgather() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  const gcore::Box & box = sys().box();

  if (sys().box().M().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");


  if (sys().numMolecules() != refSys().numMolecules())
    throw gromos::Exception("Gather problem",
          "Number of molecules in reference and frame are not the same.");
  for (int m = 0; m < sys().numMolecules(); ++m) {
    Molecule &mol = sys().mol(m);
    Molecule &refMol = refSys().mol(m);
    if (mol.numPos() != refMol.numPos()) {
      std::ostringstream msg;
      msg << "Number of positions in molecule " << m + 1 << " in reference and frame"
              " are not the same.";
      throw gromos::Exception("Gather problem", msg.str());
    }
    for (int a = 0; a < mol.numPos(); ++a) {
      // gather to the reference
      mol.pos(a) = nearestImage(refMol.pos(a), mol.pos(a), box);
      // set the current frame as the reference for the next frame
      refMol.pos(a) = mol.pos(a);
    }
  }
  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refSol = sys().sol(0);
  if (sol.numPos() != refSol.numPos()) {
    throw gromos::Exception("Gather problem", "Number of solvent positions in "
            "reference and frame are not the same.");
  }
  for (int a = 0; a < sol.numPos(); ++a) {
    // gather to the reference
    refSol.pos(a) = nearestImage(refSol.pos(a), sol.pos(a), box);
    // set the current frame as the reference for the next frame;
    refSol.pos(a) = sol.pos(a);
  }
}
