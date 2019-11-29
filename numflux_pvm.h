#ifndef NUMFLUX_PVM_H
#define NUMFLUX_PVM_H

#include <vector>

namespace global { extern int int_form; }
std::vector<std::vector<double>> numflux_pvm(std::vector<double>& ul, std::vector<double>& ur, std::vector<double>& Sl, std::vector<double>& Sr, int idx_q, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Al, std::vector<std::vector<double>>& Qurmul, std::vector<double>& Qc, std::vector<double>& S);

#endif
