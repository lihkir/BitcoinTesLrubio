#ifndef SAVE_SOLUTION_H
#define SAVE_SOLUTION_H

#include "containers.h"

void save_solution(Matrix<double> &u0, int N_rows, int N_cols, int imex_type, int convec_type, double T, int idx_q);

#endif