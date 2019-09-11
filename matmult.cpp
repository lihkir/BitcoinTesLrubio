#include "matmult.h"

void matmult(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &p)
{
    for (unsigned int i = 0; i < p.size(); i++)
        for (unsigned int j = 0; j < p[i].size(); j++)
            for (unsigned int k = 0; k < p.size(); k++)
                p[i][j] = p[i][j] + A[i][k]*x[k][j];
}