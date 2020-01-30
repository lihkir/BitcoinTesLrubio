#ifndef RKIMEX_H
#define RKIMEX_H

#include "containers.h"

namespace global { extern int idx_q; }

void rkimex(Matrix<double> &A, Vector<double> &b, Matrix<double> &Ah, Vector<double> &bh, Matrix<double> &u0, Vector<double> Ta, double cfl, double L, int convec_type0, int imex_type);
void do_rkimex(Matrix<double> &A, Vector<double> &b, Matrix<double> &Ah, Vector<double> &bh, Matrix<double> &u, double h, double dt, int maxits, double tol);
void do_lirkimex(Matrix<double> &A, Vector<double> &b, Matrix<double> &At, Matrix<double> &u, double h, double dt);

#endif
