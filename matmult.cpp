#include "matmult.h"

void matmult(Matrix<double>& A, Matrix<double>& Qurmul, Matrix<double>& fh)
{
	int N = size(Qurmul, 1); int K = size(Qurmul, 2);
	
	for (int k = 1; k <= K; k++)
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
				fh(i, k) = fh(i, k) + A(i, j)*Qurmul(j, k);
}