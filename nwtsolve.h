#ifndef NWTSOLVE_H
#define NWTSOLVE_H

#include <vector>

namespace global { extern double beta; }

void nwtsolve(double a, std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &u, double h, double dt, int maxits, double tol);

#endif
