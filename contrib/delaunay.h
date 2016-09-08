using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

class Delaunay   {
 public:
  Delaunay( const AtomSpecifier   * atoms);
  Delaunay( const AtomSpecifier   * atoms,  int i,  int j,  int k);
  Delaunay();
  ~Delaunay();
  Delaunay &operator=(const Delaunay &c);
  void set_vertices ( int i,  int j,  int k);
  void set_cog();
  void set_angle (const Boundary *  pbc , const System & sys, const Vec &  COG );
  void set_atoms (const  AtomSpecifier  * atoms);
  void set_identity (const  int  & count);
  double Calc_Angle (const Boundary *  pbc , const System & sys, const Vec &  COG,  int i);
  void writepdb(ofstream &out, int count);
  // accessors
  const std::vector< int> &  get_vertices();
  const gmath::Vec & get_cog();
  const double & get_angle();
  const int & get_identity();
 private:
  std::vector< int> vertices;
  gmath::Vec cog;
  double angle;
  int identity;
  const AtomSpecifier * Atoms;
};





