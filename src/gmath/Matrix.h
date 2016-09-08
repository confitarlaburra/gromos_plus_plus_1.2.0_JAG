// gmath_Matrix.h

#ifndef INCLUDED_GMATH_MATRIX
#define INCLUDED_GMATH_MATRIX

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_CASSERT
#include <cassert>
#define INCLUDED_CASSERT
#endif

#include <gsl/gsl_matrix.h>

namespace gmath{

  class Vec;
  /**
   * Class Matrix
   * Class that contains some basic functions for matrices
   *
   * The GROMOS++ matrix has some basic functionality. It uses the GSL
   * library for more elaburate taks
   * 
   * @class Matrix
   * @author R. Buergi, M.A. Kastenholz
   * @ingroup gmath
   */
  class Matrix{
    double **d_val;
    int d_rows, d_columns;
    
  public:
    /**
     * Matrix constructor, all elements are initialize to value
     * @param rows, columns Number of rows and columns in the matrix
     * @param value Value to which all elements are initialized
     */
    Matrix(int rows = 3, int columns = 3, double value = 0.0);
    /**
     * Matrix copy constructor
     * @param & Matrix to be copied
     */
    Matrix(const Matrix &);
    /**
     * Matrix constructor. The matrix (3x3) is constructed as the dyadic
     * product of two vectors
     * @param v,w vectors to form the matrix
     */
    Matrix(const Vec &v, const Vec &w);
    /**
     * Matrix constructor. The columns of the matrix (3x3) are formed by 
     * three vectors
     * @param u,v,w vectors to form the matrix
     */
    Matrix(const Vec &u, const Vec &v, const Vec &w);
    /**
     * Matrix deconstructor
     */    
    ~Matrix();


    // Methods
    /** 
     * Copy operator, copy matrix one into the other
     */
    Matrix &operator=(const Matrix &);
    /**
     * operator to perform a single value decomposition
     */
    Matrix luDecomp()const;
    /**
     * operator to perform an inversion using LU decomp.
     */
    Matrix invert()const;
    /**
     * an operator to diagonalise a symmetric matrix and return the 
     * eigenvalues.
     * @param eigenValues An array that is returned with the eigenvalues
     * @return The eigenvectors of the matrix
     */
    Matrix diagonaliseSymmetric(double *eigenValues, bool sort=true);
      // diagonalise a symmetric Matrix and return eigenvalues.
    /**
     * operator to calculate the determinant of a matrix
     */
    double det()const;

    /**
     * operator to calculate the determinant of a 3 by 3 matrix
     */

    double fastdet3X3Matrix()const;

    /**
     * operator to calculate the determinant of a 3 by 3 gsl-matrix
     */

    double fastdet3X3Matrix(gsl_matrix &gsl_mat);


    // operators
    /**
     * return the transpose of a matrix
     */
    Matrix transpose()const;
    /**
     * Operator that changes the sign of all the elements
     */
    Matrix operator-()const;
    /**
     * Operator to add another matrix to your matrix
     */
    Matrix &operator+=(const Matrix &mat);
    /**
     * Operator to substract another matrix from your matrix
     */
    Matrix &operator-=(const Matrix &mat);
    /**
     * Operator to multiply your matrix with a scalar
     */
    Matrix &operator*=(double d);

    // Accessors
    /**
     * Accessor that gives you the i, j element of the matrix as a const
     */
    double operator()(int i, int j)const;
    /**
     * Accessor that gives you the i, j element of the matrix
     */
    double &operator()(int i, int j);
    /**
     * Accessor that gives you the number of rows
     */
    int rows()const;
    /**
     * Accessor that gives you the number of columns
     */ 
    int columns()const;

    // Exception
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): gromos::Exception("Matrix",what){}
    };

  };

  // inline functions & free operators

  inline Matrix operator+(const Matrix &mat1, const Matrix &mat2){
    assert(mat1.rows()==mat2.rows() && mat1.columns()==mat2.columns());
    Matrix temp(mat1);
    temp+=mat2;
    return temp;
  }

  inline Matrix operator-(const Matrix &mat1, const Matrix &mat2){
    assert(mat1.rows()==mat2.rows() && mat1.columns()==mat2.columns());
    Matrix temp(mat1);
    temp-=mat2;
    return temp;
  }

  inline Matrix operator*(double d, const Matrix &m){
    Matrix temp(m);
    temp*=d;
    return temp;
  }

  Matrix operator*(const Matrix &m1, const Matrix &m2);

  Vec operator*(const Matrix &m, const Vec &v);

  inline Matrix &Matrix::operator+=(const Matrix &mat){
    assert(rows()==mat.rows() && columns()==mat.columns());
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	d_val[i][j]+=mat.d_val[i][j];
    return *this;
  }

  inline Matrix &Matrix::operator-=(const Matrix &mat){
    assert(rows()==mat.rows() && columns()==mat.columns());
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	d_val[i][j]-=mat.d_val[i][j];
    return *this;
  }

  inline Matrix &Matrix::operator*=(double d){
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	d_val[i][j]*=d;
    return *this;
  }

  inline Matrix Matrix::operator-()const{
    Matrix temp(*this);
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	temp.d_val[i][j]=-temp.d_val[i][j];
    return temp;
  }
    
  inline double Matrix::operator()(int i, int j)const{
    return d_val[i][j];
  }

  inline double &Matrix::operator()(int i, int j){
    return d_val[i][j];
  }

  inline int Matrix::rows()const{
    return d_rows;
  }
  
  inline int Matrix::columns()const{
    return d_columns;
  }

}


#endif
