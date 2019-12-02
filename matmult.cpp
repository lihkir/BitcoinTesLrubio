#include "matmult.h"
#include "utilities.h"

void matmult(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Qurmul, std::vector<std::vector<double>>& fh)
{
	for (unsigned int i = 0; i < fh.size(); i++)
		for (unsigned int j = 0; j < fh[i].size(); j++)
			fh[i][j] = fh[i][j] + dot_product(A[i], col(Qurmul, j));
}