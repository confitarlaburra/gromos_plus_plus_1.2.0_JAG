// gcore_AtomTopology.h
#ifndef INCLUDED_GCORE_ATOMTOPOLOGY
#define INCLUDED_GCORE_ATOMTOPOLOGY

#ifndef INDLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore {

  class AtomTopology_i;
  class Exclusion;

  /**
   * Class AtomTopology
   * Purpose: contains all atomic properties for an atom
   *
   * Description:
   * All topological information that can be assigned to a single atom
   * is contained in an AtomTopology
   *
   * @class AtomTopology
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::Exclusion
   * @sa gcore::MoleculeTopology
   */
  class AtomTopology {
    AtomTopology_i *d_this;

  public:
    // Constructors
    /**
     * AtomTopology constructor
     */
    AtomTopology();
    /**
     * AtomTopology copyconstructor
     * @param & AtomTopology to be copied
     */
    AtomTopology(const AtomTopology&);

    /**
     * AtomTopology deconstructor
     */
    ~AtomTopology();

    // Methods
    /**
     * Member operator = copies one AtomTopology into the other
     * @param & AtomTopology to be copied
     */
    AtomTopology & operator=(const AtomTopology&);
    /**
     * Member function to set the Integer Atom Code of the atom
     */
    void setIac(int);
    /**
     * Member function to set the ChargeGroup of the atom (0 or 1)
     */
    void setChargeGroup(int);
    /**
     * Member function to set the charge of the atom
     */
    void setCharge(double);
    /**
     * Member function to set the mass of the atom
     */
    void setMass(double);
    /**
     * Member function to set the name of the atom
     */
    void setName(const std::string&);
    /**
     * Member function to set the exclusions of the atom
     * @param e A set of exclusions of type gcore::Exclusion
     */
    void setExclusion(const Exclusion &e);
    /**
     * Member function to set the 1,4 neighbours of the atom
     * @param e A set of 1,4 neighbours of type gcore::Exclusion
     */
    void setExclusion14(const Exclusion &e);
    /**
     * Member function to set the Vanderwaals radius of the atom. This 
     * is not standard gromos96 topological information, but is used by 
     * some programs
     */
    void setradius(double);
    /**
     * Member function to set whether the atom is a hydrogen atom
     */
    void setH(bool);
    /**
     * Member function to set whether the atom is polarisable
     */
    void setPolarisable(bool);
    /**
     * Member function to set the polarisablity
     * @param a the polarisability
     */
    void setPolarisability(double a);
    /**
     * Member function to set the charge-on-spring charge
     * @param c the size of the COS charge
     */
    void setCosCharge(double c);
    /**
     * Member function to set the damping level
     * @param l the damping electric field offset
     */
    void setDampingLevel(double l);
    /**
     * Member function to set the damping power
     * @param p the damping power
     */
    void setDampingPower(double p);
    /**
     * Member function to set the gamma of the off-site polarisation centre
     * @param p the damping power
     */
    void setPoloffsiteGamma(double g);
    /**
     * Member function to set the first atom for the off-site
     * polarisation centre construction
     * @param p the damping power
     */
    void setPoloffsiteI(int i);
    /**
     * Member function to set the decond atom for the off-site
     * polarisation centre construction
     * @param p the damping power
     */
    void setPoloffsiteJ(int j);
    /**
     * Member function to set whether the atom is coarse grained
     */
    void setCoarseGrained(bool);
    /**
     * Member function to set the coarse grain factor of an atom
     */
    void setCGFactor(int);


    // Accessors
    /**
     * accessor, returns the Integer Atom Code of the atom
     */
    int iac()const;
    /**
     * accessor, returns the chargeGroup code of the atom (0 or 1)
     */
    int chargeGroup()const;
    /**
     * accessor, returns the charge of the atom
     */
    double charge()const;
    /**
     * accessor, returns the mass of the atom
     */
    double mass()const;
    /** 
     * accessor, returns the name of the atom
     */
    const std::string &name()const;
    /** 
     * accessor, returns the exclusions of the atom (const)
     */
    const Exclusion &exclusion()const;
    /** 
     * accessor, returns the 1,4 neighbours of the atom (const)
     */
    const Exclusion &exclusion14()const;
    /** 
     * accessor, returns the exclusions of the atom
     */
    Exclusion &exclusion();
    /** 
     * accessor, returns the 1,4 neighbours of the atom
     */
    Exclusion &exclusion14();
    /**
     * accessor, returns the Vanderwaals radius of the atom
     */
    double radius()const;
    /**
     * accessor, returns whether the atom is a hydrogen atom. 
     * Can also be used to flag your atoms for other things
     * @sa MoleculeTopology setHmass setHtype
     */
    bool isH()const;
    /**
     * accessor, returns whether the atom is polarisable
     */
    bool isPolarisable()const;
    /**
     * accessor, returns the polarisability of the atom
     */
    double polarisability()const;
    /**
     * accessor, returns the size of the COS charge connected to the atom
     */
    double cosCharge()const;
    /**
     * accessor, returns the damping level electric field offset of the 
     * atom
     */
    double dampingLevel()const;
    /**
     * accessor, returns the damping power
     */
    double dampingPower()const;
    /**
     * accessor, returns the gamma of the off-site polarisation
     */
    double poloffsiteGamma()const;
    /**
     * accessor, returns the first atom defining the off-site polarisation
     * centre
     */
    int poloffsiteI()const;
    /**
     * accessor, returns the decond atom defining the off-site polarisation
     * centre
     */
    int poloffsiteJ()const;
    /**
     * accessor, returns whether the atom is coarse grained
     */
    bool isCoarseGrained()const;
    /**
     * accessor, returns the coarse grain factor
     */
    int cg_factor()const;
  };
}
#endif

