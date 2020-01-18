#ifndef CONVEC_GLF_H
#define CONVEC_GLF_H

#include "containers.h"

void convec_glf(Matrix<double> &v, double h, Matrix<double> &K);
void numflux_glf(Matrix<double> &v, Matrix<double> &fh);

#endif
