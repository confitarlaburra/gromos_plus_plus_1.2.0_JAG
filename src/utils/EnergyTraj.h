// utils_CheckTopo.h

// Class that runs some basic checks on a molecule topology
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_MAP
#include <map>
#define INCLUDED_MAP
#endif
#ifndef INCLUDED_GMATH_EXPRESSION
#include "../gmath/Expression.h"
#define INCLUDED_GMATH_EXPRESSION
#endif
#ifndef INCLUDED_GIO_GINSTREAM
#include "../gio/Ginstream.h"
#define INCLUDED_GIO_GINSTREAM
#endif

namespace utils{
  /**
   * Class EnergyIndex
   * A class that is an index to a data array
   */
  class EnergyIndex
  {
  public:
    /**
     * block index, to find the data
     */
    int block;
    /**
     * index i to the element of the block
     */
    int i;
    /**
     * index j to the element of the block
     */
    int j;
    /**
     * a list of the dependencies for an expression
     */
    std::vector<EnergyIndex> dep;
    /**
     * an operator to see if two EnergyIndex are the same
     */
    bool operator==(EnergyIndex const& ei)const;
    
  };
  
  /**
   * Class EnergyBlock
   * The definition a data block type
   */
  class EnergyBlock 
  {
  public:
    /**
     * enum to determine the type of the dimensionality
     */
    enum type_enum{fixed, var, size, matsize, block};
    /**
     * type for the first dimension
     */
    type_enum first_type;
    /**
     * type for the second dimension
     */
    type_enum second_type;
    /**
     * size (or index to the size) of the first dimension
     */
    int first_size;
    /**
     * size (or index to the size) of the second dimension
     */
    int second_size;
    /**
     * index to the datablock where the data is actually stored
     */
    int blockindex;
    /**
     * name of the block
     */
    std::string blockname;
    /**
     * function to read the block from an instream
     */
    int read(gio::Ginstream &gin, std::istringstream & iss, 
	 std::vector<std::vector<std::vector<double > > > & blocks,
	 std::vector<int> & sizes);
    /**
     * constructor
     */
    EnergyBlock(std::string s, int b, type_enum t1, int i1,
		type_enum t2=fixed, int i2=1)
    {
      blockname=s;
      blockindex=b;
      first_type=t1;
      second_type=t2;
      first_size=i1;
      second_size=i2;
      /**
      if(first_type==size || first_type==matsize){
	std::cout << "size " << i1;
      }
      else {
	std::cout << "block (" << blockindex << "):";
	if(first_type==fixed)
	  std::cout << i1;
	else
	  std::cout << "var " << i1;
	std::cout << " , ";
	if(second_type==fixed)
	  std::cout << i2;
	else
	  std::cout << "var " << i2;
      }
      std::cout << std::endl;
      */
    }
    EnergyBlock(std::string s)
    {
      first_type=block;
      second_type=fixed;
      blockname=s;
    }
    

  };
  
  
    
  
  /**
   * Class EnergyTraj
   * A class that contains all energy terms from (free) energy trajectories
   * 
   * The EnergyTraj class stores the energy terms which can be accessed via 
   * (user defined) names. It also allows the user to define how properties 
   * should be calculated from other values, defined through a 
   * gmath::Expression 
   *
   * @class EnergyTraj
   * @author B.C. Oostenbrink
   * @ingroup utils
   */
  class EnergyTraj
    {
      /**
       * A vector of doubles that will contain all energy trajectory data
       *
       */
      std::vector<std::vector<std::vector<double> > > d_data;
      /**
       * A vector of EnergyIndex
       */
      std::vector<utils::EnergyIndex> d_index;
      /**
       * A vector of vectors of EnergyBlock (for every file a vector)
       */
      std::vector<std::vector<utils::EnergyBlock> > d_blocks;
      /**
       * A vector of sizes
       */
      std::vector<int> d_size;
      /**
       * A vector of expressions that define how a property should be computed
       *
       * The properties that need to be computed get internally a negative 
       * index.
       */
      std::vector<gmath::Expression> d_e;
      /**
       * A vector that stores the result of the properties that need to be
       * calculated so that a second access does not recalculate.
       */
      std::vector<double> d_calc;
      /**
       * A vector of booleans that keeps track of which properties have
       * been calculated already and which one need to be recalculated.
       */
      std::vector<bool> d_recalc;
      /**
       * A type definition for the map
       */
      typedef std::map<std::string, EnergyIndex>::value_type MP;
      /** 
       * A map that links the name of a property to its index.
       * Since we do not know the number of energy groups in advance, these
       * energies are not put in the map, but are handled seperately in the 
       * index function
       */
      std::map<std::string, EnergyIndex> d_index_map;
      /**
       * A map that links the name of a size variable to its size index
       */
      std::map<std::string, int> d_size_map;
      /**
       * A map that links the block name to the block index
       */
      std::map<std::string, int> d_block_map;
      /**
       * A map that links a file type name to an index
       */
      std::map<std::string, int> d_file_map;
      
      /**
       * a counter for the number of energy frames that have been read.
       */
      int d_en_frame;
      /**
       * a counter for the number of free energy frames that have been read.
       */
      int d_fr_frame;
      /**
       * A constant that will be returned if one tries to access
       * an unknown element
       */
      int unknownvariable;
      /**
       * A string to keep track of the first kind of block in the energy 
       * trajectory
       */
      std::string en_first;
      /**
       * A string to keep track of the first kind of block in the free
       * energy trajectory
       */
      std::string fr_first;
      
    public:
      /**
       * Constructor
       *
       * This constructor constructs a standard Energy Trajectory wich has
       * all the elements that are defined on page III-56 of the GROMOS96 
       * manual.
       */
      EnergyTraj();
      /**
       * Accessor to get the element in the data set that is referred to with
       * the name s
       *
       * The name of an element can be the name according to the (free) energy
       * files (page III-56 of the GROMOS manual) or a user defined name (see
       * function addKnown
       */
      double operator[](std::string s);
      /**
       * function to read in one frame from the energy trajectory.
       */
      int read_frame(gio::Ginstream& is, std::string file_type);
      /**
       * function to teach the class about a new block of data
       */
      void addBlock(std::string s, std::string file_type);
      /**
       * function to teach the class about a new property
       *
       * A new property can be either a direct mapping of an existing known 
       * property, or it can be an expression, that can be calculated from 
       * previously known properties.
       * @param string s The name of the new property
       * @param string v An expression (or name) of the existing properties
       *                 that it refers to.
       */
      void addKnown(std::string s, std::string v);
      /**
       * function to teach the class about a new constant
       */
      void addConstant(std::string s, double f);
      
      /**
       * A function that returns the index-number of the internal data array
       * belonging to the property with name s
       */
      EnergyIndex index(std::string s);
      /**
       * A function that writes out the definitions and connections of all 
       * known properties. All known properties are listed in alphabetical 
       * order, followed by their definition. This can be just a mapping to the
       * data array, or it can be a user defined property for which the
       * definition and the dependencies are also listed.
       */
      void write_map(std::ostream& os = std::cout);
    private:
      /**
       * A function that gives the value of a property as function of the
       * index
       */
      double value(EnergyIndex const & ei);
      
      /**
       * A function that backtracks the first name that can be found,
       * which refers to the index number i
       */
      std::string back_index(EnergyIndex const &i);
      /**
       * A function to Tokenize a string into a vector of strings
       */
      void Tokenize(const std::string& str,
		    std::vector<std::string>& tokens,
		    const std::string& delimiters = " ");
      /**
       * A function to initialize the free energy trajectory. Set some
       * variables and learn about the standard names;
       */
      void init();
    };
}



      
  
