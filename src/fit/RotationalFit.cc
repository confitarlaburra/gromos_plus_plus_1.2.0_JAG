// fit_RotationalFit.cc
//includes explicit calls to gsl now

#include <cassert>
#include "RotationalFit.h"
#include "Reference.h"
#include "PositionUtils.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Matrix.h"
#include "../gmath/Vec.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <iostream>

using fit::RotationalFit;
using fit::RotationalFit_i;
using fit::Reference;
using gcore::System;
using gcore::Molecule;
using gmath::Matrix;
using gmath::Vec;

// constructs the rotation Matrix
static void rotationMatrix(gmath::Matrix *mat, const gcore::System &mol, const fit::Reference &w);

RotationalFit::RotationalFit(Reference *w){
  d_ref=w;
  PositionUtils::shiftToCog(&w->sys(),*w);
}

RotationalFit::~RotationalFit(){}

void RotationalFit::fit(gcore::System *sys)const{
  //check that ref has enough atoms
  int count=0;
  for(int m =0; m<sys->numMolecules(); ++m)
    for(int a=0; a < sys->mol(m).numAtoms(); ++a)
      if(d_ref->weight(m,a)!=0.0)
	count++;
  if(count<=3)
    throw RotationalFit::Exception("At least 4 atoms required for rotational fit\n");
  
	  
  PositionUtils::shiftToCog(sys,*d_ref);
  Matrix rot(3,3);
  rotationMatrix(&rot,*sys,*d_ref);
  PositionUtils::rotate(sys,rot);
}

static void rotationMatrix(Matrix *mat, const System &sys, const Reference &r){

  const System &ref = r.sys();
  
  
  Matrix U(3,3,0);
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      for(int m=0;m<ref.numMolecules();++m)
  	for(int n=0;n<ref.mol(m).numAtoms();++n)
  	  if(r.weight(m,n))
	    U(i,j)+=r.weight(m,n)*sys.mol(m).pos(n)[i]*ref.mol(m).pos(n)[j];

  double det=U.fastdet3X3Matrix();

  int signU = ( det>0 ? 1 : -1);
  
  
  gsl_matrix * omega = gsl_matrix_alloc (6, 6);
  gsl_matrix_set_zero (omega);
  
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      gsl_matrix_set (omega, i, j+3,U(i,j));
      gsl_matrix_set (omega, i+3, j,U(j,i));
    }
  }
  
  
  double *eigenvals = new double [6];
  
  gsl_vector *eval = gsl_vector_alloc (6);
  gsl_matrix *evec = gsl_matrix_alloc (6,6);
  
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (6);

  gsl_eigen_symmv (omega, eval, evec, w);
  
  gsl_eigen_symmv_free(w);

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  
  Matrix Omega(6,6,0);
  for (int i=0; i < 6; ++i){
    eigenvals[i] = gsl_vector_get(eval, i);
    for (int j=0; j < 6; ++j){
      Omega(i,j)=gsl_matrix_get(evec, i, j);
    }
  }
  
  gsl_matrix_free (omega);
  gsl_matrix_free (evec);
  gsl_vector_free (eval);

  if(det<0 && fabs(eigenvals[1] - eigenvals[2]) < 1.0e-8){

    std::cerr << "determinant = " << det << "\n"
	      << "eigenval[0] = " << eigenvals[0] << "\n"
	      << "eigenval[1] = " << eigenvals[1] << "\n"
	      << "eigenval[2] = " << eigenvals[2] << "\n" << std::endl;
    
    throw RotationalFit::Exception("Rotation matrix degenerate!");
  }
  
  // Extract vectors from Omega.  
  Omega *= sqrt(2.0);
  Vec k1(Omega(0,0), Omega(1,0), Omega(2,0));
  Vec k2(Omega(0,1), Omega(1,1), Omega(2,1));
  Vec k3(Omega(0,2), Omega(1,2), Omega(2,2));
  Vec h1(Omega(3,0), Omega(4,0), Omega(5,0));
  Vec h2(Omega(3,1), Omega(4,1), Omega(5,1));
  Vec h3(Omega(3,2), Omega(4,2), Omega(5,2));

  double spat = h1.dot(h2.cross(h3));
  
  // turn 3rd vectors
  if(spat<0){
    h3=-h3;
    k3=-k3;
  }

  *mat = Matrix(h1,k1) + Matrix(h2,k2) + signU*Matrix(h3,k3);

  delete[] eigenvals;
 
}
