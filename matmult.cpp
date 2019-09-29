#include "matmult.h"
#include "utilities.h"

void matmult(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& p)
{
	for (unsigned int i = 0; i < p.size(); i++)
		for (unsigned int j = 0; j < p[i].size(); j++)
			p[i][j] = p[i][j] + DotProduct(A[i], Col(x, j));
}