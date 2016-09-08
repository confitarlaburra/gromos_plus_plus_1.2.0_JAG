#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "EnergyTraj.h"
#include "../gmath/Expression.h"
#include "../gromos/Exception.h"
#include "../gio/Ginstream.h"

#include <cassert>

using utils::EnergyTraj;
using utils::EnergyBlock;
using utils::EnergyIndex;


EnergyTraj::EnergyTraj()
{
  init();
}

void EnergyTraj::init()
{
  unknownvariable=-1000;

  d_en_frame=0;
  d_fr_frame=0;

  d_block_map["CONSTANTS"]=0;
  d_data.resize(1);
  
}

utils::EnergyIndex EnergyTraj::index(std::string s)
{  
  // first slow attempt
  std::string::size_type it1=s.find("["), it2;
  std::string name=s;

  utils::EnergyIndex ei;
  
  if(it1 != std::string::npos){
    name=s.substr(0,it1);
  }
  //  std::cout << "accessing name " << name << std::endl;
  
  if(d_block_map.count(name)){
    ei.block=d_block_map[name];
    
    it2=s.find("]");
    ei.i = atoi(s.substr(it1+1,it2).c_str())-1;
    ei.j = 0;
    it1=s.find("[",it1+1);
    it2=s.find("]",it2+1);
    if(it1!=std::string::npos && it2!=std::string::npos){
      ei.j = atoi(s.substr(it1+1, it2).c_str())-1;
    }
  }
  else{
    if(d_index_map.count(name)){
      ei = d_index_map[name];
    }
    else
      ei.i = unknownvariable;
  }
  
  //  std::cout << "\twhich is " << ei.block<< " : " << ei.i << " : " << ei.j << std::endl;
  
  return ei;
}
double EnergyTraj::value(EnergyIndex const & ei)
{
  if(ei.block>=0){
    if(ei.block > int(d_data.size())  || 
       ei.i > int(d_data[ei.block].size())  ||
       ei.j > int(d_data[ei.block][ei.i].size())){
      
      
      throw gromos::Exception("EnergyTraj", 
			      "Trying to access an unknown variable ");
      return 0.0;
    }
    else
      return d_data[ei.block][ei.i][ei.j];
  }
  else{
    
    // it is an expression that we know how to handle
    int ind = -ei.block -1;
    //    std::cout << "we calculate an expression " << ind << std::endl;
    //    d_e[ind].writeExpression(std::cout);
    //    d_e[ind].writeExpressionValue(std::cout);

    if(ind > int(d_recalc.size()))
      throw gromos::Exception("EnergyTraj", 
			      "Trying to calculate an expression that is "
			      "unknown to me");
    
    if(d_recalc[ind]){
      std::vector<double> v;
      for(unsigned int k=0; k<ei.dep.size(); k++)
	v.push_back(value(ei.dep[k]));
      d_e[ind].setValues(v);   
      //      d_e[ind].writeExpressionValue(std::cout);

      d_recalc[ind]=false;
      d_calc[ind]=d_e[ind].value();
    }
    return d_calc[ind];
  }
}

double EnergyTraj::operator[](std::string s)
{
  utils::EnergyIndex ei = index(s);
  if(ei.i == unknownvariable) 
    throw gromos::Exception("EnergyTracetory", 
			    "Trying to access an unknown variable " +s);
  return value(ei);
}

void EnergyTraj::addConstant(std::string s, double f)
{
  //  std::cout << "adding constant " << s << " " << f << std::endl;
  
  //check if it exists already
  EnergyIndex ei=index(s);
  //  std::cout << "ei.i " << ei.i << std::endl;
  
  if(ei.i!=unknownvariable){
    //    std::cout << "setting known constant" << std::endl;
    
    d_data[ei.block][ei.i][ei.j]=f;
  }
  else{
    //    std::cout << "setting new constant " << std::endl;
    
    ei.block = d_block_map["CONSTANTS"];
    //    std::cout << ei.block << std::endl;

    ei.i = d_data[ei.block].size();
    ei.j = 0;
    //    std::cout << ei.i << " " << ei.j << std::endl;
    d_data[ei.block].resize(ei.i+1);
    d_data[ei.block][ei.i].resize(ei.j+1);
    
    d_data[ei.block][ei.i][ei.j]=f;
    d_index_map[s]=ei;
  }
}

void EnergyTraj::addKnown(std::string s, std::string e)
{
  //  std::cout << "entering addKnown " << s << " : " << e << std::endl;

  // Tokenize the expression string
  std::vector<std::string> v;
  Tokenize(e,v);
  //  std::cout << v.size() << std::endl;
  
  // if it is just a one-to-one mapping, find it and add it to the map;
  if(v.size()==1){
    EnergyIndex ei = index(v[0]);
    //    std::cout << ei.block << " " << ei.i << " " << ei.j << std::endl;
    
    if(ei.i != unknownvariable){
      d_index_map[s]=ei;
    }
    else{
      //it must be a constant
      double f;
      std::stringstream is(v[0]);
      is >> f;
      addConstant(s,f);
    }
  }
  else{
    // we have an expression
    std::ostringstream os;
    int varcount = 1;
    EnergyIndex ei=index(s);
    ei.dep.resize(0);
    
    for(unsigned int j=0; j< v.size(); j++){
      EnergyIndex var=index(v[j]);
      if(var.i!=unknownvariable){
	os << "a" << varcount;
	ei.dep.push_back(var);
	varcount++;
      }
      else os << v[j];
      os << " ";
    }

    if(ei.i==unknownvariable){
      gmath::Expression e(os.str());

      d_e.push_back(e);
      d_calc.push_back(0.0);
      d_recalc.push_back(true);
      ei.block=-1*int(d_e.size());
      ei.i=0;
    }
    else{
      int ind=-1*ei.block-1;
      
      d_e[ind].setExpression(os.str());
      d_calc[ind]=0;
      d_recalc[ind]=true;
      
    }
    d_index_map[s]=ei;
  }
  
  
  /*
  }
  else{

    }
    
    gmath::Expression e(os.str());
    
    if(!known){
      d_e.push_back(e);
      d_dep.push_back(dep);
      d_calc.push_back(0.0);
      d_recalc.push_back(true);
      d_map.insert(MP(s, -1*d_e.size()));
    }
    else{
      d_e[i].setExpression(os.str());
      d_dep[i].resize(0);
      for(unsigned int k=0; k<dep.size(); k++)
	d_dep[i].push_back(dep[k]);
    }
  }
  */
}

int EnergyTraj::read_frame(gio::Ginstream& gin, std::string file_type)
{
  int fileindex;
  if(d_file_map.count(file_type))
    fileindex=d_file_map[file_type];
  else
    throw gromos::Exception("EnergyTraj", "Don't know how to read from file"
			    " type " + file_type);
  //std::cout << "fileindex " << fileindex << std::endl;

  std::istringstream is;
  
  //std::cout << "blocksize " << d_blocks[fileindex].size() << std::endl;
  
  for(unsigned int i=0; i< d_blocks[fileindex].size(); i++){
    if(d_blocks[fileindex][i].read(gin, is, d_data, d_size))
      return 1;
  }
  // set everything that needs to be calculated to be recalculated
  for(unsigned int i=0; i<d_recalc.size(); i++) d_recalc[i]=true;

  return 0;
}

int EnergyBlock::read(gio::Ginstream &gin, std::istringstream &iss, 
		       std::vector<std::vector<std::vector<double > > > & data,
		       std::vector<int> & size)
{
  std::string s;
  std::vector<std::string> buffer;
  //std::cout << "block name " << blockname << std::endl;
  
  if(first_type==block){
    gin.getblock(buffer);
    if(gin.stream().eof()) return 1;
    
    if(buffer[0]!=blockname)
      throw gromos::Exception("EnergyTraj", "Block " + blockname + " expected. "
			      "Got " + buffer[0]);

    //std::cout << "reading block " << blockname << std::endl;
    
    gio::concatenate(buffer.begin()+1, buffer.end()-1, s);
    iss.clear();
    iss.str(s);
  }
  else if(first_type==EnergyBlock::size){
    int i=0;
    if(!( iss >> i) )
      throw gromos::Exception("EnergyTraj", "Tried to read an integer for " +
			      blockname);
    size[first_size]=i;
    size[second_size]=i*(i+1)/2;
  }  
  else {
    unsigned int si=first_size;
    unsigned int sj=second_size;
    if(first_type==var)  si = size[first_size];
    if(second_type==var) sj = size[second_size];
    //std::cout << "resizing " << si << " " << sj << std::endl;
    if(data[blockindex].size() < si)
      data[blockindex].resize(si);
    for(unsigned int i=0; i<si; i++) 
      if(data[blockindex][i].size() < sj) 
	data[blockindex][i].resize(sj);
    
    for(unsigned int i=0; i<si; i++)
      for(unsigned int j=0; j<sj; j++)
	if(!(iss >> data[blockindex][i][j]))
	  throw gromos::Exception("EnergyTraj", "Tried to read a value in " +
				 blockname);
    
  }
  return 0;
}

bool EnergyIndex::operator==(EnergyIndex const& ei)const
{
  return block==ei.block && i==ei.i && j==ei.j;
}

std::string EnergyTraj::back_index(utils::EnergyIndex const& ei)
{
  // first search the blocks	
  // in the case of a constant or an expression, we would rather have our 
  // own definition and not a sizetype with the same index, so we don't search 
  // for them in the blocks
  if(ei.block>0){
    
    for(unsigned int i=0; i<d_blocks.size(); i++){
      for(unsigned int j=0; j<d_blocks[i].size(); j++){
	if(d_blocks[i][j].blockindex==ei.block){
	  std::ostringstream os;
	  os << d_blocks[i][j].blockname;
	  os << "[" << ei.i+1 << "]";
	  if(d_blocks[i][j].second_type!=EnergyBlock::fixed ||
	     d_blocks[i][j].second_size!=1) 
	    os << "[" << ei.j+1 << "]";
	  return os.str();
	}
      }
    }
  }
  
  // then search the map
  std::map<std::string, EnergyIndex>::const_iterator iter=d_index_map.begin(), 
    to=d_index_map.end();
  for(; iter!=to; ++iter)
    if(iter->second == ei) return iter->first;

  return "unknown";
}

void EnergyTraj::write_map(std::ostream& os)
{
  std::map<std::string, EnergyIndex>::const_iterator iter=d_index_map.begin(), 
    to=d_index_map.end();
  for(;iter!=to; ++iter){
    os << std::setw(12) << iter->first << " = ";
    os << "data[" << iter->second.block << "][" << iter->second.i 
       << "][" << iter->second.j << "]";
    std::string nm=back_index(iter->second);
    if(iter->first != nm)
      os << "  (" << nm << ")"; 
    else if(iter->second.block==0)
      os << "  (constant: " << value(iter->second) << ")";
    
    os << std::endl;
    if(iter->second.block < 0){
      os << std::setw(15) << "= ";
      int i=-1-iter->second.block;
      d_e[i].writeExpression(os);
      if(iter->second.dep.size()!=0){
	os << std::setw(20) << "with:" << std::endl;
	for(unsigned int j=0; j<iter->second.dep.size(); j++)
	  os << std::setw(20) << "a" << j+1 << " = data[" 
	     << iter->second.dep[j].block << "]["
	     << iter->second.dep[j].i     << "]["
	     << iter->second.dep[j].j     << "] (" 
	     << this->back_index(iter->second.dep[j])
	     << ")" << std::endl;
      }
    }
  }
}

void EnergyTraj::Tokenize(const std::string& str,
			  std::vector<std::string>& tokens,
			  const std::string& delimiters)
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

void EnergyTraj::addBlock(std::string s, std::string file_type)
{
  int fileindex=0;

  // std::cerr << d_file_map.size() << std::endl;
  // std::cerr << "filetype = " << file_type << std::endl;
  
  if(d_file_map.count(file_type)){
    fileindex=d_file_map[file_type];
    // std::cerr << "file_type => fileindex = " << fileindex << std::endl;
  }
  else{
    // std::cerr << d_file_map.size() << std::endl;
  
    const int d_file_map_size = d_file_map.size();
    
    d_file_map[file_type]=d_file_map_size;

    fileindex=d_file_map[file_type];
    d_blocks.resize(d_file_map.size());
  }
  // std::cerr  << "fileindex " << fileindex << std::endl;
  // std::cerr << "d_blocks: "<< d_blocks.size() << std::endl;
  
  std::vector<std::string> t;
  Tokenize(s, t);
  std::transform(t[0].begin(), t[0].end(), t[0].begin(), tolower);
  int i, j;
  EnergyBlock::type_enum i_type=EnergyBlock::fixed;
  EnergyBlock::type_enum j_type=EnergyBlock::fixed;
  
  if(t[0]=="subblock"){
    // std::cerr << "subblock" << std::endl;
    
    std::istringstream is(t[2]);
    if(t.size()!=4)
      throw gromos::Exception("EnergyTraj", "Not enough block parameters");
    
    if(!(is >> i)) {
      // std::cerr << "i not a number" << std::endl;

      if(d_size_map.count(t[2]))
	i=d_size_map[t[2]];
      else if(d_size_map.count("matrix_"+t[2]))
	i=d_size_map["matrix_"+t[2]];
      else
	throw gromos::Exception("EnergyTraj", "Variable name "+t[2]+
				  " in block size not defined");
      i_type=EnergyBlock::var;
    }

    is.clear();
    is.str(t[3]);
    if(!(is >> j)) {
      // std::cerr << "j not a number" << std::endl;

      if(d_size_map.count(t[3]))
	j=d_size_map[t[3]];
      else if(d_size_map.count("matrix_"+t[3]))
	j=d_size_map["matrix_"+t[3]];
      else
	throw gromos::Exception("EnergyTraj", "Variable name "+t[3]+
				" in block size not defined");
      
      j_type=EnergyBlock::var;
    }

    if(!d_block_map.count(t[1])){
      const int d_block_map_size = d_block_map.size();
      d_block_map[t[1]]=d_block_map_size;
      d_data.resize(d_block_map.size());
    }
    
    assert(fileindex >=0 && fileindex < int(d_blocks.size()));
    d_blocks[fileindex].push_back(EnergyBlock(t[1], d_block_map[t[1]], 
					      i_type, i, j_type, j));
    
  }
  else if(t[0] == "size"){
    // std::cerr << "size" << std::endl;
    
    if(t.size() != 2)
      throw gromos::Exception("EnergyTraj", "Not enough size parameters");
    if(!d_size_map.count(t[1])){

      const int d_size_map_size = d_size_map.size();
      d_size_map[t[1]]=d_size_map_size;

      //std::cout << "sizemap " << d_size_map.size() << std::endl;
      // also store the associated matsize
      const int d_matrix_size_map_size = d_size_map.size();
      d_size_map["matrix_"+t[1]]=d_matrix_size_map_size;

      d_size.resize(d_size_map.size());
    }
    
    // std::cerr << "size " << d_size.size() << std::endl;
    // std::cerr << "fileindex " << fileindex << std::endl;
    
    assert(fileindex >=0 && fileindex < int(d_blocks.size()));
    d_blocks[fileindex].push_back(EnergyBlock(t[1], 0, 
					      EnergyBlock::size, 
					      d_size_map[t[1]],
					      EnergyBlock::matsize,
					      d_size_map["matrix_"+t[1]]));
  }
  else if(t[0] == "block"){
    // std::cerr << "block" << std::endl;
    if(t.size() != 2)
      throw gromos::Exception("EnergyTraj", "Not enough block parameters");

    assert(fileindex >=0 && fileindex < int(d_blocks.size()));
    d_blocks[fileindex].push_back(EnergyBlock(t[1]));
  }
  else
    throw gromos::Exception("EnergyTraj", "Don't know how to handle keyword "+
			    t[0]);
}


    

