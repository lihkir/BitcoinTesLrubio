#ifndef RKIMEX_H
#define RKIMEX_H

#include <vector>

void rkimex(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& Ah, std::vector<double>& bh, std::vector<std::vector<double>>& u0, std::vector<double> Ta, double cfl, double L, int convec_type, int imex_type, int test);
void do_rkimex(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& Ah, std::vector<double>& bh, std::vector<std::vector<double>>& u, double h, double dt, int maxits, double tol);
double do_lirkimex(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& At, std::vector<std::vector<double>>& u, double h, double dt);

#endif
