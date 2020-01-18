#ifndef HORNERM_H
#define HORNERM_H

#include <stdio.h>
#include "matmult.h"
#include "matmult_dec.h"
#include "containers.h"

namespace global { extern int int_form; }

void hornerm(Matrix<double>& A, Matrix<double>& dif, Matrix<double>& Qc, Matrix<double>& fh, Matrix<double>& Qurmul);

#endif 
