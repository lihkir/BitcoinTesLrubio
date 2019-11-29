#ifndef NEWTON_MATRIX_H
#define NEWTON_MATRIX_H

#include <vector>

void newton_matrix(std::vector<std::vector<double>> &u, double h, double a, std::map<int, std::vector<std::vector<double>>> &D, std::map<int, std::vector<std::vector<double>>> &L, std::map<int, std::vector<std::vector<double>>> &U);

#endif
