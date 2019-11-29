#ifndef RESIDUAL_H
#define RESIDUAL_H

#include <vector>
#include <map>

void residual(std::map<int, std::vector<std::vector<double>>>& D, std::map<int, std::vector<std::vector<double>>>& L, std::map<int, std::vector<std::vector<double>>>& U, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& r);

#endif
