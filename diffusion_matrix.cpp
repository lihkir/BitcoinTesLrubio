#include <math.h>
#include "diffusion_tensor.h"
#include "diffusion_matrix.h"

void diffusion_matrix(Matrix<double>& u, double h, double a, Block<double> &D, Block<double> &L, Block<double> &U)
{
	int n = size(u, 2);
	int m = size(u, 1);

	double z;
	double f = -a*0.5/(h*h);

	Block<double> *ptB = new Block<double>(n, m, m); Block<double> &B = *ptB;
	diffusion_tensor(u, B);

	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			for (int p = 1; p <= m; p++)
			{
				D(j, p, i) = 0;
            	L(j, p, i) = 0;
            	U(j, p, i) = 0;
			}
		}
	}

	for (int p = 1; p <= m; p++)
	{
		for (int j = 1; j <= n; j++)
		{
			if (j > 1)
			{
				for (int i = 1; i <= m; i++)
				{
					/** K(p, j) = K(p, j)-(B(p, i, j)+B(p, i, j-1))*(u(i, j)-u(i, j-1)) **/
                	z = B(j, p, i) + B(j - 1, p, i);
					D(j, p, i) = D(j, p, i) - f*z;
                	L(j - 1, p, i) = L(j - 1, p, i) + f*z;
				}
			}
        	if (j <  n)
			{
				for (int i = 1; i <= m; i++)
				{
					/** K(p, j)=K(p, j)+(B(p, i, j+1)+B(p, i, j))*(u(i, j+1)-u(i, j)) **/
                	z = B(j + 1, p, i) + B(j, p, i);
                	U(j, p, i) = U(j, p, i) + f*z;
                	D(j, p, i) = D(j, p, i) - f*z;
				}
			}
		}
	}

	for (int j = 1; j <= n; j++)
		for (int i = 1; i <= m; i++)
			D(j, i, i) = D(j, i, i) + 1;

	delete ptB;
}