// gcore_CrossDihedral.cc

#include "CrossDihedral.h"
#include "Dihedral.h"
#include <new>

using gcore::CrossDihedral;
using gcore::Dihedral;

CrossDihedral::CrossDihedral(int a, int b, int c, int d, int e, int f, int g, int h) {
  Dihedral d1(a, b, c, d), d2(e, f, g, h);
  if (d1 < d2) {
    d_a[0] = d1[0];
    d_a[1] = d1[1];
    d_a[2] = d1[2];
    d_a[3] = d1[3];
    d_a[4] = d2[0];
    d_a[5] = d2[1];
    d_a[6] = d2[2];
    d_a[7] = d2[3];
  } else {
    d_a[0] = d2[0];
    d_a[1] = d2[1];
    d_a[2] = d2[2];
    d_a[3] = d2[3];
    d_a[4] = d1[0];
    d_a[5] = d1[1];
    d_a[6] = d1[2];
    d_a[7] = d1[3];
  }
  d_type = -1;
}

CrossDihedral::CrossDihedral(const CrossDihedral &a) {
  for (unsigned int i = 0; i < 8; ++i)
    d_a[i] = a.d_a[i];
  d_type = a.d_type;
}

CrossDihedral &CrossDihedral::operator=(const CrossDihedral &b) {
  if (this != &b) {
    this->CrossDihedral::~CrossDihedral();
    new(this) CrossDihedral(b);
  }
  return *this;
}

int gcore::operator<(const CrossDihedral &a, const CrossDihedral &b) {
  const Dihedral a1(a[0], a[1], a[2], a[3]), b1(b[0], b[1], b[2], b[3]);
  if (a1 < b1)
    return 1;
  if (b1 < a1)
    return 0;
  // a1 == b1
  const Dihedral a2(a[4], a[5], a[6], a[7]), b2(b[4], b[5], b[6], b[7]);
  if (a2 < b2)
    return 1;
  if (b2 < a2)
    return 0;
  // a2 == b2
  return a.type() < b.type();
}

