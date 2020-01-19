#include "matmult_dec.h"
#include "update_utilities.h"

void matmult_dec(Matrix<double>& A, Matrix<double>& Qurmul, Matrix<double>& fh)
{
	int N = Qurmul.rows();
	int K = Qurmul.cols();

    Matrix<double> *ptQurmulk = new Matrix<double>(N, 1); Matrix<double> &Qurmulk = *ptQurmulk;

    if (N == 1)
        matmult(A, Qurmul, fh);
    else
    {
	    for (int k = 1; k <= K; k++)
        {
            /** p=p+diag(Adec(:,1))*x **/
            for (int i = 1; i <= N; i++) 
                fh(i, k) = fh(i, k) + A(i, 1)*Qurmul(i, k);        
            /** s=sum(x); p=p+Adec(:,2)*s **/
            get_col(Qurmulk, Qurmul, k);
            double s = sum1D(Qurmulk);
            for (int i = 1; i <= N; i++) 
                fh(i, k) = fh(i, k) + A(i, 2)*s;     
        }        
    }
    delete ptQurmulk;
}