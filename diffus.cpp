#include "diffus.h"
#include "diffusion_tensor.h"
#include "test_cases.h"
#include "update_utilities.h"

void diffus(Matrix<double>& u1, double h, Matrix<double> &K)
{
	struct test_cases* pt_test = get_tests();
	
	int n = size(u1, 2) - 2*pt_test->gc;
	int m = size(u1, 1);

	Matrix<double> *ptu = new Matrix<double>(m, n); Matrix<double> &u = *ptu;
	update_inside(u, u1, pt_test->gc);

	Block<double> *ptB = new Block<double>(n, m, m); Block<double> &B = *ptB;
	diffusion_tensor(u, B);

	for (int p = 1; p <= m; p++)
	{
		for (int j = 1; j <= n; j++)
		{
			K(p, j) = 0;
			if (j > 1)
			{
				for (int i = 1; i <= m; i++)
					K(p, j) = K(p, j) - (B(j, p, i) + B(j - 1, p, i))*(u(i, j) - u(i, j-1));
			}
			if (j < n)
			{
				for (int i = 1; i <= m; i++)
					K(p, j) = K(p, j) + (B(j + 1, p, i) + B(j, p, i))*(u(i, j + 1) - u(i, j));
			}
            K(p, j) = 0.5*K(p, j)/(h*h);
		}
	}
	delete ptB;
	delete ptu;
}