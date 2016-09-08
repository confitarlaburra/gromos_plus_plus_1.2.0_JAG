#include "delaunay.h"

/*
using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
*/

Delaunay::Delaunay(const AtomSpecifier  * atoms) : 
  vertices(3,0),
  cog(0.0,0.0,0.0),
  angle(0.0)
{
  Atoms = atoms;
}

Delaunay::Delaunay(const AtomSpecifier   * atoms,  int i,  int j,  int k) :
  cog(0.0,0.0,0.0),
  angle(0.0)
{
  Atoms= atoms;
  vertices[0]=i;
  vertices[1]=j;
  vertices[2]=k;
}

Delaunay::Delaunay() {}
Delaunay::~Delaunay() {}


Delaunay &Delaunay::operator=(const Delaunay &c){
  if (this != &c) {
    vertices = c.vertices;
    cog= c.cog;
    angle= c.angle;
    Atoms=c.Atoms;
  }
  return *this;
}
void Delaunay::set_vertices( int i,  int j,  int k) {
  vertices[0]=i;
  vertices[1]=j;
  vertices[2]=k;
}

void Delaunay::set_cog(){
  cog=(Atoms->pos(vertices[0])+ Atoms->pos(vertices[1])+ Atoms->pos(vertices[2]))*(1.0/3.0);
}

void Delaunay::set_angle(const Boundary *  pbc , const System & sys, const Vec &  COG ) {
  
  //Angle between cog,COG,vertices[0]
  gmath::Vec tmpA, tmpB;
  tmpA = cog  - pbc->nearestImage(cog,COG, sys.box());
  tmpB = Atoms->pos(vertices[0])- pbc->nearestImage(Atoms->pos(vertices[0]), COG, sys.box());
  angle = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;
  
}

void Delaunay::set_identity (const  int  & count) {
  identity=count;
}


// Computes angle beween any vertix and the cog of the triangle
double Delaunay::Calc_Angle (const Boundary *  pbc , const System & sys, const Vec &  COG,  int i){
  //Angle between cog,COG,atom i
  gmath::Vec tmpA, tmpB;
  double angle_vertix;
  tmpA = cog  - pbc->nearestImage(cog,COG, sys.box());
  tmpB = Atoms->pos(i)- pbc->nearestImage(Atoms->pos(i), COG, sys.box());
  angle = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;
  return angle_vertix;
}

void Delaunay::set_atoms ( const AtomSpecifier  * atoms) {
  Atoms= atoms;
}

// writes the vertices in pdb format with CONECT
void Delaunay::writepdb(ofstream &out, int count) {

   out.setf(ios::fixed, ios::floatfield);
   out.setf(ios::unitbuf);
   out.precision(3);
   vector<int>::const_iterator iter=vertices.begin();
   vector<int>::const_iterator to=vertices.end();
   for (; iter != to; iter++) {
     out << "HETATM";
     out.setf(ios::right, ios::adjustfield);
     out << setw(5) << count;
     out.setf(ios::left, ios::adjustfield);
     out << "  " <<setw(4) << (Atoms->name(*iter)).c_str();
     out << setw(4) << (Atoms->resname(*iter)).c_str() <<" "<< setw(4) << (Atoms->resnum(*iter)) << "    ";
     out.setf(ios::right, ios::adjustfield);
     
     out  << setw(8) << ((Atoms->pos(*iter))[0])*10
	  << setw(8) << ((Atoms->pos(*iter))[1])*10
	  << setw(8) << ((Atoms->pos(*iter))[2])*10
	  << "  1.00  0.00" << endl;
     count++;
   }

    out << "CONECT "
	<< setw(4) << count - 2 << setw(10) << count -1 << endl;

    out << "CONECT "
	<< setw(4) << count - 1 << setw(10) << count  << endl;
}

//accesors
const std::vector<int> & Delaunay::get_vertices(){
  return vertices;
};
const gmath::Vec &  Delaunay::get_cog(){
  return cog;
};
const double &  Delaunay::get_angle(){
  return angle;
};

const int &  Delaunay::get_identity(){
  return identity;
};




