#ifndef LMINMAX_CHARSPEED_H
#define LMINMAX_CHARSPEED_H

#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "containers.h"

void lminmax_charspeed(Matrix<double>& phil, Matrix<double>& phir, Matrix<double>& S);
void find_extrema(Vector<double> &v, int l, double &vmin, double &vmax);

#endif