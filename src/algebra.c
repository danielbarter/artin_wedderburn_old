#include "algebra.h"

double get_multiplication_3(Algebra *algebra, int i, int j, int k) {
  int d = algebra->dimension;
  int index = k + d * j + d * d * i;
  return algebra->multiplication[k];
}

void set_multiplication_3(Algebra *algebra, double v, int i, int j, int k) {
  int d = algebra->dimension;
  int index = k + d * j + d * d * i;
  algebra->multiplication[k] = v;
}

double *get_multiplication_2(Algebra *algebra, int i, int j) {
    int d = algebra->dimension;
    int offset = d * j + d * d * i;
    return algebra->multiplication + offset;
}

void set_multiplication_2(Algebra *algebra, double *vp, int i, int j) {
    int d = algebra->dimension;
    int offset = d * j + d * d * i;
    double *base_pointer = algebra->multiplication + offset;

    int k; // loop index
    for (k = 0; k < d; k++) {
      base_pointer[k] = vp[k];
    }
}

double *get_multiplication_1(Algebra *algebra, int i) {
    int d = algebra->dimension;
    int offset = d * d * i;
    return algebra->multiplication + offset;
}

void set_multiplication_1(Algebra *algebra, double *vp, int i) {
    int d = algebra->dimension;
    int offset = d * d * i;
    double *base_pointer = algebra->multiplication + offset;

    int jk; // loop index
    for (jk = 0; jk < d * d; jk++) {
      base_pointer[jk] = vp[jk];
    }
}
