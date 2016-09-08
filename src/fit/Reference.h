// fit_Reference.h

#ifndef INCLUDED_FIT_WEIGHTS
#define INCLUDED_FIT_WEIGHTS

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{
  class System;
}

namespace utils{
  class AtomSpecifier;
}

namespace fit{

  class Reference_i;
  /**
   * Class Reference
   * The reference class keeps weights for the atoms we are interested in
   *
   * This is used for e.g. a rotational fit or an rmsd calculation
   *
   * @class Reference
   * @author R. Buergi
   * @ingroup fit
   * @sa fit::RotationalFit
   * @sa fit::TranslationalFit
   * @sa utils::Rmsd
   */
  class Reference{
    Reference_i *d_this;

    // not implemented
    Reference();
    Reference(const Reference &);
    Reference &operator=(const Reference &);

  public:
    /**
     * Reference constructor
     * A reference always needs information about the system
     */
    Reference(gcore::System *sys);
    /**
     * Reference deconstructor
     */
    ~Reference();
    
    // methods
    /**
     * Method to add a certain class of atoms to the reference
     * All atoms that have a certain name get the same weight (1)
     * @param mol The number of the molecule you are interested in
     * @param name The name of the atoms that should be added
     */
    void addClass(int mol, const std::string &name);
    /**
     * Method to add a specific atom to the Reference
     * @param m The number of the Molecule the Atom belongs to
     * @param i The number of the Atom in this Molecule
     */
    void addAtom(int m, int i);    
    /**
     * Method to add an AtomSpecifier to the Reference
     * @param as the AtomSpecifier to be added
     */
    void addAtomSpecifier(utils::AtomSpecifier as);
    /**
     * Method to give equal weight to all non-zero elements
     */
    void rescale();
    /**
     * Method to set the weight of a specific atom to a certain value
     * @param m The number of the molecule of this atom
     * @param i The number of the atom
     * @param w The weight that the atom should get
     */
    void setWeight(int m, int i, double w);
    /**
     * method to normalize all weights, to have a total summed up value of
     * 1
     */
    void normalise();
    /**
     * A method that makes an integer list containing all atom numbers that
     * have a specific name
     * @param molecule The number of the molecule you are interested in
     * @param atomtype The name of the atoms that you want to have
     * @param poslist is returned with a vector of integers containing the 
     *                numbers of all atoms with this name
     */
   void makePosList (const gcore::System &,int molecule, const std::string &atomtype, std::vector<int> &poslist);

    // accessor
   /**
    * Accessor that returns the system on which the Reference is based
    */
    gcore::System &sys();
    /**
     * Accessor that returns the system on which the Reference is based 
     * as a const
     */
    const gcore::System &sys()const;
    // return reference to system.
    /**
     * Accessor that returns the weight of atom i in molecule m
     */
    double weight(int m, int i)const;
  };
}

#endif
