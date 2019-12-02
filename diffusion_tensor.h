#ifndef DIFFUSION_TENSOR_H
#define DIFFUSION_TENSOR_H

#include <vector>
#include <map>

std::map<int, std::vector<std::vector<double>>> diffusion_tensor(std::vector<std::vector<double>>& u);

#endif