#ifndef NUMFLUX_PVM_H
#define NUMFLUX_PVM_H

#include "containers.h"

namespace global { extern int int_form; }

void numflux_pvm(Matrix<double> &ul, Matrix<double> &ur, Matrix<double> &Sl, Matrix<double> &Sr, int idx_q, Matrix<double> &fh, Matrix<double> &A, Matrix<double> &Al, Matrix<double> &fl, Matrix<double> &fr, Matrix<double> &ua, Matrix<double> &dif, Matrix<double> &Qurmul, Matrix<double> &Qc, Matrix<double> &S);

#endif
