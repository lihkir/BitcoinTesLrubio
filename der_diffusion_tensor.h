#ifndef DER_DIFFUSION_TENSOR_H
#define DER_DIFFUSION_TENSOR_H

#include <vector>
#include <map>

namespace global { extern double D0; }

std::map<int, std::vector<std::vector<double>>> der_diffusion_tensor(std::vector<std::vector<double>>& u);

#endif