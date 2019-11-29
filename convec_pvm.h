#ifndef CONVEC_PVM_H
#define CONVEC_PVM_H

#include <vector>

namespace global { extern int idx_q; }

void convec_pvm(std::vector<std::vector<double>> &v, double h, std::vector<std::vector<double>>&K);

#endif
