#ifndef CONVEC_GLF_H
#define CONVEC_GLF_H

#include <vector>

void convec_glf(std::vector<std::vector<double>> &v, double h, std::vector<std::vector<double>> &K);
void numflux_glf(std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &fh);

#endif
