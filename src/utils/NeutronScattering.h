namespace utils {

  class iNS;

  /**
   * Class NS (Neutron Scattering) calculates the scattering intensities ...
   *
   * And some more descriptions here please...
   *
   * @class NS
   * @ingroup utils
   * @author A. Eichenberger
   */
  class NS {

  private:

    /**
     * A pointer to the implementation class containing the data
     */
    class iNS *d_this;
    /**
     * The default constructor: should not be used since it is likely
     * to forget the setting of the ponters and iterators by hand later on...
     */
    NS(void);

  public:

    /**
     * Constructor
     */
    NS(gcore::System *sys, const args::Arguments *args);

    /**
     * Destructor to delete the data (implementation class)
     */
    ~NS(void);

    /**
     * Sets the number of grid points of the resulting spectrum.
     * @param grid the number of grid points
     */
    void setGrid(int grid);
    /**
     * Sets the cutoff used in the RDF calculations
     * @param cut the number of grid points
     */
    void setCut(double cut);
    /**
     * Sets the maximum Q-value of the resulting spectrum.
     * @param Qmax the maximum Q-value
     */
    void setQmax(int Qmax);
    /**
     * Sets the atoms to be considered as centre atoms
     */
    int addAtoms(std::string s);
    /**
     * Gets all combination of centre to with atom types and resets the lengths
     * of the corresponding vectors to that length
     */
    int getCombinations(void);
    /**
     * Gets/sets the weighting (pre)factors for the inter- and intra-molecular
     * partial structure factors. Make sure the centre and with atoms are
     * set before using this function
     */
    void getWeights(void);
    /**
     * Checks the lengths and initialisation state of the members of NS to
     * test if it is ready for calculation of the scattering intensities
     */
    void check(void);
    /**
     * Sets a system variable to all dependant variables
     * @param sys a gromos system (gcore::System)
     */
    void setSystem(gcore::System *sys);
    /**
     * Set the centre and with atoms of all d_rdf members
     */
    void setRDFatoms();
    /**
     * Calculates all inter-molecular g(r) which are needed to build up the
     * structure factors S(Q) as well as the intensities I(Q)
     */
    void calcRDFsInterAll();
    /**
     * Reads the scattering lengths from file
     */
    void readScattlen(std::string fname);
    /**
     * Reads the atom pair positional root-mean-square position deviation
     */
    void readSigma(std::string fname);
    /**
     * Prints the NEUTRONSCATTERING block
     */
    void print(std::ostream &os);
    /**
     * The sinc function
     */
    double sinc(double x);
    /**
     * Calculate the intra-molecular structure factors
     */
    void calcSintra(void);
    /**
     * Calculate the elastic intra-molecular structure factors
     */
    void calcSintraElastic(void);
    /**
     * Calculate the inter-molecular structure factors
     */
    void calcSinter(void);
    /**
     * Calculates the scattering intensity I(Q) based on the partial
     * iter- and intra-molecular strucutre factors
     */
    void calcIntensity(void);
    /**
     * Prints the neutron scattering intensity
     */
    void printIntensity(std::ostream &os);
    /**
     * Prints the partial strucutre factors
     */
    void printS(std::ostream &os);
    /**
     * prints the different radial distribution functions
     */
    void printRDFs(std::ostream &os);
    /**
     * Prints the intra-molecular partial structure factors
     */
     void printSintra(std::ostream &os);
     /**
     * Prints the inter-molecular partial structure factors
     */
     void printSinter(std::ostream &os);


  };

} /* end of namespace utils */
