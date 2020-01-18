#ifndef NEWTON_MATRIX_H
#define NEWTON_MATRIX_H

#include "containers.h"

void newton_matrix(Matrix<double> &u, double h, double a, Block<double> &D, Block<double> &L, Block<double> &U);

#endif
