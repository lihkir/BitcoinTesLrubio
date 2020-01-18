#ifndef UPDATE_UTILITIES_H
#define UPDATE_UTILITIES_H

#include "containers.h"

void update_col(const Matrix<double>& u, const Matrix<double>& c, int j);
void copy_block(Block<double>& u, Block<double>& v);
void update_inside(Matrix<double>& u, Matrix<double>& v, int gc);
void loop_over(Matrix<double> &u, Matrix<double> &v, int N_col, int gcl, int gcr);
void get_col(Matrix<double>& c, Matrix<double>& u, int j);

#endif
