#ifndef DIFFUSION_MATRIX_H
#define DIFFUSION_MATRIX_H

#include <vector>
#include <map>

void diffusion_matrix(std::vector<std::vector<double>>& u, double h, double a, std::map<int, std::vector<std::vector<double>>> &D, std::map<int, std::vector<std::vector<double>>> &L, std::map<int, std::vector<std::vector<double>>> &U);

#endif
