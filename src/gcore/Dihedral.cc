// gcore_Dihedral.cc

#include "Dihedral.h"
#include <new>
#include <iostream>

using gcore::Dihedral;

Dihedral::Dihedral(int a, int b, int c, int d, bool warn){
  if(b<c){
    d_a[0]=a; 
    d_a[1]=b; 
    d_a[2]=c; 
    d_a[3]=d;
  } 
  else{
    d_a[0]=d;
    d_a[1]=c;
    d_a[2]=b;
    d_a[3]=a;
    if (warn) {
      std::cerr << "NOTE: order of atoms changed in dihedral:\n";
      std::cerr << "      " << a + 1 << "," << b + 1 << "," << c + 1 << "," << d + 1 << " -> "
              << d + 1 << "," << c + 1 << "," << b + 1 << "," << a + 1 << std::endl;
    }
  }
  d_type=-1;
}

Dihedral::Dihedral(const Dihedral &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_a[2]=a.d_a[2];
  d_a[3]=a.d_a[3];
  d_type=a.d_type;
}

Dihedral &Dihedral::operator=(const Dihedral &b){
  if(this != &b){
    this->Dihedral::~Dihedral();
    new(this) Dihedral(b);
  }
  return *this;
}

int gcore::operator<(const Dihedral &a, const Dihedral &b){
  return (a[1]<b[1]||(a[1]==b[1]&&a[2]<b[2])
	  ||(a[1]==b[1]&&a[2]==b[2]&&a[0]<b[0])
	  ||(a[1]==b[1]&&a[2]==b[2]&&a[0]==b[0]&&a[3]<b[3])
	  ||(a[1]==b[1]&&a[2]==b[2]&&a[0]==b[0]&&a[3]==b[3]&&a.type()<b.type()));
}
