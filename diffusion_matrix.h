#ifndef DIFFUSION_MATRIX_H
#define DIFFUSION_MATRIX_H

#include "containers.h"

void diffusion_matrix(Matrix<double>& u, double h, double a, Block<double> &D, Block<double> &L, Block<double> &U);

#endif
