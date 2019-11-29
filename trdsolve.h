#ifndef TRDSOLVE_H
#define TRDSOLVE_H

#include <vector>
#include <map>

void trdsolve(std::map<int, std::vector<std::vector<double>>>& A1, std::map<int, std::vector<std::vector<double>>>& B1, std::map<int, std::vector<std::vector<double>>>& C1, std::vector<std::vector<double>>& b);
void lutrb(std::map<int, std::vector<std::vector<double>>> &A, std::map<int, std::vector<std::vector<double>>> &B, std::map<int, std::vector<std::vector<double>>> &C);
void nplu(std::vector<std::vector<double>> &A);
void fwdsolve(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
void bwdsolve(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
void utrsolve(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
void maxpy(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &y);

#endif
