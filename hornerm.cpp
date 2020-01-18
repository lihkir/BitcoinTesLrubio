#include "hornerm.h"

void hornerm(Matrix<double>& A, Matrix<double>& dif, Matrix<double>& Qc, Matrix<double>& fh, Matrix<double>& Qurmul)
{
	int N = size(dif, 1);
	int K = size(dif, 2);
	int n = size(Qc, 1); 

	for (int k = 1; k <= K; k++)
		for (int i = 1; i <= N; i++)
			Qurmul(i, k) = Qc(n)*dif(i, k);

	for (int k = 1; k <= K; k++)
	{
		for (int j = n - 1; j >= (1); j-=1)
		{
			if (Qc(j) != 0)
				for (int i = 1; i <= N; i++) fh(i, k) = Qc(j)*dif(i, k);
			else
				for (int i = 1; i <= N; i++) fh(i, k)=0;
			
			/** p=p+A*pAx = a*x+A*pAx **/
			if (global::int_form != 5) {
				matmult(A, Qurmul, fh);
			} else {
				matmult_dec(A, Qurmul, fh);
			}
			/** pAx= p = a*x+A*pAx **/
			for (int i = 1; i <= N; i++) Qurmul(i, k) = fh(i, k);
		}
	}
}