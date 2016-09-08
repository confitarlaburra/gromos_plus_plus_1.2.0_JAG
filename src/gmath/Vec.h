// gmath_Vec

#ifndef INCLUDED_GMATH_VEC
#define INCLUDED_GMATH_VEC

#ifndef INCLUDED_CASSERT
#include <cassert>
#define INCLUDED_CASSERT
#endif

#ifndef INCLUDED_CMATH
#include <cmath>
#define INCLUDED_CMATH
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gmath{
  /**
   * Class Vec
   * gromos++ definition of a 3 element vector
   *
   * The gromos vector class has some specific properties for 
   * vectors with length 3
   *
   * @class Vec
   * @author R. Buergi
   * @ingroup gmath
   */ 
  class Vec{
    double d_v[3];

  public:
    /**
     * Vec constructor
     * @param i, j, k the elements that form the vector
     */
    Vec(double i=0, double j=0, double k=0);
    /**
     * Vec copy constructor
     */
    Vec(const Vec &v);
    /**
     * Vec deconstructor
     */
    ~Vec(){}
    // Methods
    /**
     * Operator to assign one vector to another
     */
    Vec &operator=(const Vec &v);
    /**
     * Operator to change the sign of all elements in a vector
     */
    Vec operator-()const;
    /**
     * Operator to add one vector to another vector
     */
    Vec &operator+=(const Vec &v);
    /**
     * Operator to substract another vector from your vector
     */
    Vec &operator-=(const Vec &v);
    /**
     * Operator to multiply the elements of two vectors with each other
     */
    Vec &operator*=(const Vec &v);
    /**
     * Operator to multiply the elements of the vecor with a scalar
     */
    Vec &operator*=(double d);
    /**
     * Operator to divide the elements of the vector by a scalar
     */
    Vec &operator/=(double d);
    /**
     * Operator to determine if one vector equals another
     */
    bool operator==(const Vec &v)const;

    /**
     * Cross (outer) product of two gmath::Vectors
     */
    Vec cross(const Vec &v)const;
    /**
     * Dot (inner) product of two gmath::Vectors
     */
    double dot(const Vec &v)const;
    /** 
     * Method that returns the normalized gmath::Vector, i.e. the Vector 
     * divided by its absolute value
     */
    Vec normalize()const;
    /**
     * The dot product of a gmath::Vector with itself (norm)
     */
    double abs2()const;
    /**
     * The sqrt of the norm of the gmath::Vector (absolute value)
     */
    double abs()const;

    
    /**
     * Accessor that returns the i-th element of a vector as a const
     */
    double operator[](int i)const;
    /**
     * Accessor that returns the i-th element of a vector
     */
    double &operator[](int i);
  };


  std::string v2s(gmath::Vec const &v);


  // Class Methods
  
  inline Vec::Vec(double i, double j, double k){
    d_v[0]=i;
    d_v[1]=j;
    d_v[2]=k;
  }
 
  inline Vec::Vec(const Vec &v){
    d_v[0]=v.d_v[0];
    d_v[1]=v.d_v[1];
    d_v[2]=v.d_v[2];
  }

  inline Vec Vec::operator-()const{
    return Vec(-d_v[0], -d_v[1], -d_v[2]);
  }
 
  inline Vec &Vec::operator=(const Vec &v){
    if(this != &v){
      d_v[0]=v.d_v[0];
      d_v[1]=v.d_v[1];
      d_v[2]=v.d_v[2];
    }   
    return *this;
  }

  inline Vec &Vec::operator+=(const Vec &v){
    d_v[0]+=v.d_v[0];
    d_v[1]+=v.d_v[1];
    d_v[2]+=v.d_v[2];
    return *this;
  }

  inline Vec &Vec::operator*=(const Vec &v){
    d_v[0]*=v.d_v[0];
    d_v[1]*=v.d_v[1];
    d_v[2]*=v.d_v[2];
    return *this;
  }

  
  inline Vec &Vec::operator-=(const Vec &v){
    d_v[0]-=v.d_v[0];
    d_v[1]-=v.d_v[1];
    d_v[2]-=v.d_v[2];
    return *this;
  }

  inline Vec &Vec::operator*=(double d){
    d_v[0]*=d;
    d_v[1]*=d;
    d_v[2]*=d;
    return *this;
  } 

  inline Vec &Vec::operator/=(double d){
    d_v[0]/=d;
    d_v[1]/=d;
    d_v[2]/=d;
    return *this;
  } 

  inline Vec Vec::cross(const Vec &v) const{
    return Vec(d_v[1]*v.d_v[2] - d_v[2]*v.d_v[1],
	       d_v[2]*v.d_v[0] - d_v[0]*v.d_v[2],
	       d_v[0]*v.d_v[1] - d_v[1]*v.d_v[0]);
  }



  inline double Vec::dot(const Vec &v) const{
    return d_v[0]*v.d_v[0] + d_v[1]*v.d_v[1] + d_v[2]*v.d_v[2];
  }

  inline double Vec::abs2()const{
    return this->dot(*this);
  }

  inline double Vec::abs()const{
    return std::sqrt(this->abs2());
  }

  inline Vec Vec::normalize()const{
    double d = std::sqrt(this->abs2());
    return Vec(d_v[0]/d, d_v[1]/d, d_v[2]/d);
  }

  inline double Vec::operator[](int i)const{
    assert(i<3);
    return d_v[i];
  }

  inline double &Vec::operator[](int i){
    assert(i<3);
    return d_v[i];
  }

  inline bool Vec::operator==(const Vec &v)const {
    return (
	    operator[](0) == v[0] &&
	    operator[](1) == v[1] &&
	    operator[](2) == v[2]
	    );
  }

  inline Vec operator+(const Vec &a, const Vec &b){
    Vec v(a);
    v+=b;
    return v;
  }

  inline Vec operator-(const Vec &a, const Vec &b){
    Vec v(a);
    v-=b;
    return v;
  }

  inline Vec operator*(const Vec &a, const Vec &b){
    Vec v(a);
    v*=b;
    return v;
  }

  inline Vec operator*(double d, const Vec &b){
    Vec v(b);
    v*=d;
    return v;
  }

  inline Vec operator*(const Vec &a, double d){
    Vec v(a);
    v*=d;
    return v;
  }

  inline Vec operator/(const Vec &a, double d){
    Vec v(a);
    v/=d;
    return v;
  }
  
}
#endif





