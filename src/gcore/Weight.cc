#include <map>
#include <string>
#include "../args/Arguments.h"
#include "../gio/Ginstream.h"
#include "../gromos/Exception.h"
#include "Weight.h"

static const double DEFAULT_WEIGHT = 1.0;

void gcore::Weight::read(const args::Arguments& args, const std::string argname) {
  args::Arguments::const_iterator it = args.lower_bound(argname),
          to = args.upper_bound(argname);
  if (it == to) {
    has_weights = false;
    return;
  }
  has_weights = true;
  for(; it != to; ++it) {
    gio::Ginstream istr(it->second);
    istr.readTitle();
    std::vector<std::string> block;
    if (block[0] != "CONFIGURATIONWEIGHTS") {
      std::ostringstream os;
      os << "File " << it->second << " did not contain a CONFIGURATIONWEIGHTS block.";
      throw gromos::Exception("Weight", os.str());
    }
    std::string block_content;
    gio::concatenate(block.begin()+1, block.end()-1, block_content);
    std::istringstream iss(block_content);
    double w;
    while(iss >> w) {
      weights.push(w);
    }
  }
}

double gcore::Weight::getWeight() const {
  if (has_weights)
    return weights.front();
  return DEFAULT_WEIGHT;
}

void gcore::Weight::pop() {
  // do not pop the first time
  if (has_weights) {
    if (!first) {
      if (weights.empty()) {
        throw gromos::Exception("Weight", "There were not enough weights provided. "
                "Please give a weight number for every configuration.");
      }
      weights.pop();
    } else {
      first = false;
    }
    sum_weights += weights.front();
  } else {
    sum_weights += DEFAULT_WEIGHT;
  }
}
