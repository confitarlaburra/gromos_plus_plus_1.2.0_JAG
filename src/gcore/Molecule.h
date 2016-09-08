// gcore_Molecule.h
/**
 * Class Molecule
 * Addition: velocity configuration added to Molecule definition;
 * Author  : gee          
 */

#ifndef INCLUDED_GCORE_MOLECULE
#define INCLUDED_GCORE_MOLECULE

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gmath{
class Vec;
}

using gmath::Vec;

/**
 * Class Molecule
 * Purpose: Contains coordinates and topology of a Molecule
 *
 * Description:
 * The gromos++ Molecule is a specific entity. It is defined as a set of 
 * atoms that are connected by bonds. The class Molecule contains both the
 * coordinates of the molecule and the topological information
 *
 * @class Molecule
 * @ingroup gcore
 * @author R. Buergi
 * @ingroup gcore
 * @sa gcore::MoleculeTopology
 * @sa gcore::System
 */
namespace gcore{

  class MoleculeTopology;

  class Molecule{
    MoleculeTopology *d_mt;
    std::vector<Vec*> d_pos;
    std::vector<Vec*> d_vel; 
    std::vector<Vec*> d_cosDisplacement;
    
    // not implemented
    Molecule();
    Molecule &operator=(const Molecule &);

  public:
    /**
     * Molecule constructor
     * @param mt is a MoleculeTopology, a new molecule (with no coordinates)
     *           will be constructed based on this MoleculeTopology
     */
    Molecule(const MoleculeTopology &mt);
    /**
     * Molecule copy constructor
     * @param & the molecule to be copied
     */
    Molecule(const Molecule &);
    /**
     * Molecule deconstructor
     */
    ~Molecule();
    
    //Accessors
    /**
     * Accessor, returns the number of atoms in the Molecule
     */
    int numAtoms()const;
    /**
     * Accessor, returns the MoleculeTopology
     */
    MoleculeTopology &topology();
    /**
     * Accessor, returns the MoleculeTopology of the Molecule as a const
     */
    const MoleculeTopology &topology() const; 
    /**
     * Accessor, returns a pointer to a vector with the coordinates 
     * of the i-th atom in the Molecule
     * @return A gmath::Vec of three coordinates
     */
    Vec &pos(int i);
     /**
     * Accessor, returns the coordinates of the i-th atom in the Molecule 
     * as a const
     * @return A gmath::Vec of three coordinates
     */
    const Vec &pos(int i)const;
    /*
     * Accessor, returns a pointer to a vector with the velocity
     * of the i-th atom in the Molecule
     * @return A gmath::Vec of three coordinates
     */
    Vec &vel(int i);

    /**
     * Accessor, returns the velocity of the i-th atom in the Molecule 
     * as a const
     * @return A gmath::Vec of three coordinates
     */
    const Vec &vel(int i)const;
    
    /*
     * Accessor, returns the charge-on-spring displacement of the i-th atom 
     * in the Molecule
     * @return A gmath::Vec of three coordinates
     */
    Vec &cosDisplacement(int i);
    /**
     * Accessor, returns the charge-on-spring displacement of the i-th atom
     * int the Molecule as a const
     * @return A gmath::Vec of three corrdinates
     */
    const Vec &cosDisplacement(int i)const;

    /**
     * function to resize the coordinate array to the number of atoms
     * in the topology
     */
    void initPos();

    /**
     * function to resize the velocity array to the number of atoms
     * in the topology
     */
    void initVel();
    
    /**
     * function to resize the COS displacement array to the number of atoms
     * in the topology
     */
    void initCosDisplacements();

    /**
     * Accessor, returns the number of position coordinates for the Molecule
     */
    int numPos()const;
    /**
     * Accessor, returns the number of velocity co-ordinates for the Molecule
     */
    int numVel()const;
    /**
     * Accesssor, returns the number of charge-on-spring displacements for the
     * Molecule
     */
    int numCosDisplacements()const;
    
  }; /* class Molecule */

  inline Vec &Molecule::pos(int i){
    assert(i < this->numPos());
    return *d_pos[i];
  }

  inline const Vec &Molecule::pos(int i)const{
    assert (i < this->numPos());
    return *d_pos[i];
  }

  inline Vec &Molecule::vel(int i){
    assert(i < this->numVel());
    return *d_vel[i];
  }
  
  inline Vec &Molecule::cosDisplacement(int i){
    assert(i < this->numCosDisplacements());
    return *d_cosDisplacement[i];
  }

  inline const Vec &Molecule::cosDisplacement(int i)const{
    assert (i < this->numCosDisplacements());
    return *d_cosDisplacement[i];
  }

  inline const Vec &Molecule::vel(int i)const{
    assert (i < this->numVel());
    return *d_vel[i];
  }
  inline int Molecule::numPos()const{
    return d_pos.size();
  }
  inline int Molecule::numVel()const{
    return d_vel.size();
  }
  inline int Molecule::numCosDisplacements()const{
    return d_cosDisplacement.size();
  }

  
} /* Namespace */ 
#endif

