#ifndef CONVEC_PVM_H
#define CONVEC_PVM_H

#include "containers.h"

namespace global { extern int idx_q; }

void convec_pvm(Matrix<double> &v, double h, Matrix<double>&K);

#endif
