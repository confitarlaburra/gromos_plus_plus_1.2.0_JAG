// args_GatherParser.h

#ifndef INCLUDED_ARGS_GATHERPARSER
#define INCLUDED_ARGS_GATHERPARSER

#include "../bound/Boundary.h"

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif


namespace gcore{
  class System;
}

namespace bound{
  class Boundary;
}

namespace args{
  class Arguments;


/**
 * Class GatherParser
 * Purpose: Parse gathering methods from args.
 * 
 *
 * Description:
 * This class is used to parse gathering methods from args.
 * By default (when no arguments are given) it will use the gather-gromos 
 * gathering method. By default it will use args["pbc"].
 *
 * 
 * @class GatherParser
 * @version $Date: Mon Jul 15 14:17:41 MEST 2002
 * @author  M.A. Kastenholz
 * @ingroup args
 * @sa args::Arguments
 * @sa args::BoundaryParser
 */
  class GatherParser{
    
    // not implemented
    GatherParser(const GatherParser&);
    GatherParser &operator=(const GatherParser&);

  public:
/**
 * GatherParser constructor.
 * Details.
 */
    GatherParser(){}
/**
 * GatherParser destructor.
 * Details.
 */
    ~GatherParser(){}
/** 
 * Constructs the class and returns a member pointer to a gathering method.
 * Method parse parses the input from args.
 * @param args Arguments from the input line.
 * @param str name of the argument string (default "pbc")
 * @return bound::Boundary::MemPtr Member pointer to the gathering method.
 * Details.
 */
    static bound::Boundary::MemPtr parse(gcore::System &sys,gcore::System &refSys,const Arguments &gathargs, const std::string &str = "pbc");

  };

}

#endif
