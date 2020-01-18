#ifndef NWTSOLVE_H
#define NWTSOLVE_H

#include "containers.h"

namespace global { extern double beta; }

void nwtsolve(double a, Matrix<double> &b, Matrix<double> &u, double h, double dt, int maxits, double tol);

#endif
