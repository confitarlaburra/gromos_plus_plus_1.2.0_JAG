// Matrix.t.cc

#include "Matrix.h"
#include "Vec.h"
#include <iostream>
#include <gsl/gsl_matrix.h>

using namespace gmath;
using namespace std;

ostream & operator<<(ostream &os, const Matrix &mat) {
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.columns(); ++j)
      cout << mat(i, j) << ' ';
    cout << endl;
  }
  return os;
}

int main() {
  Matrix mat(3, 3);

  mat(0, 0) = 0;
  mat(0, 1) = 1;
  mat(0, 2) = 1;
  mat(1, 0) = 1;
  mat(1, 1) = 0;
  mat(1, 2) = 1;
  mat(2, 0) = 1;
  mat(2, 1) = 1;
  mat(2, 2) = 0;

  cout << "Original Matrix:\n" << mat;

  Matrix cmat(mat);
  cout << "Copied Matrix:\n" << cmat;

  int *index;
  index = new (int [mat.rows()]);

  cout << "LU Decomosition: Matrix\n";
  Matrix lumatt(mat);
  Matrix lumat = lumatt.luDecomp();
  cout << lumat;
  cout << lumat << endl;

  Matrix U(lumat);
  U(1, 0) = 0;
  U(2, 0) = 0;
  U(2, 1) = 0;

  Matrix L(lumat - U);
  L(0, 0) = 1;
  L(1, 1) = 1;
  L(2, 2) = 1;

  lumat = L*U;

  cout << "L*U:\n" << lumat;

  cout << "Determinant: " << mat.det() << endl;

  cout << "Determinant 3 by 3: " << mat.fastdet3X3Matrix() << endl;

  gsl_matrix * omega = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_zero(omega);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      gsl_matrix_set(omega, i, j, mat(i, j));
    }
  }

  cout << "Determinant 3 by 3 gsl: " << mat.fastdet3X3Matrix() << endl;

  gsl_matrix_free(omega);

  lumat = mat;
  cout << "After =:\n" << lumat;
  double eigen[3];
  Matrix ei = lumat.diagonaliseSymmetric(eigen);
  cout << "Diagonalised matrix: \n" << ei
          << "Eigenvalues: " << eigen[0] << ' '
          << eigen[1] << ' ' << eigen[2] << endl;

  Vec v(1, 1, 1);
  cout << "Product with (1,1,1):\n";
  Vec w = mat*v;
  cout << w[0] << ' ' << w[1] << ' ' << w[2] << endl;

  return 0;
}
