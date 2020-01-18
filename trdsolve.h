#ifndef TRDSOLVE_H
#define TRDSOLVE_H

#include "containers.h"

void trdsolve(Block<double> &A1, Block<double> &B1, Block<double> &C1, Matrix<double> &b);
void lutrb(Block<double> &A, Block<double> &B, Block<double> &C);
void nplu(Matrix<double> &A);
void fwdsolve(Matrix<double> &A, Matrix<double> &B);
void bwdsolve(Matrix<double> &A, Matrix<double> &B);
void utrsolve(Matrix<double> &A, Matrix<double> &B);
void maxpy(Matrix<double> &A, Matrix<double> &x, Matrix<double> &y);

#endif
