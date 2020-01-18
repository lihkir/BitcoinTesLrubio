#include "apply_diffus.h"
#include "test_cases.h"
#include "diffusion_tensor.h"
#include "update_utilities.h"

void apply_diffus(Matrix<double> &ut, Matrix<double> &uh, double h, Matrix<double> &K)
{
	struct test_cases* pt_test = get_tests();
	
	int n = size(uh, 2) - 2*pt_test->gc;
	int m = size(uh, 1);

	Matrix<double> u(m, n);
	Matrix<double> ut0(m, n);

	update_inside(u, uh, pt_test->gc);
	update_inside(ut0, ut, pt_test->gc);
	
	Block<double>  B(n, m, m);
	diffusion_tensor(ut0, B);

	for (int p = 1; p <= m; p++)
	{
		for (int j = 1; j <= n; j++)
		{
			K(p, j) = 0;
			if (j > 1)
			{
				for (int i = 1; i <= m; i++)
					K(p, j)=K(p, j)-(B(j, p, i) + B(j - 1, p, i))*(u(i, j)-u(i, j-1));
			}
			if (j <  n)
			{
				for (int i = 1; i <= m; i++)
					K(p, j) = K(p, j)+(B(j + 1, p, i)+B(j, p, i))*(u(i, j+1)-u(i, j));
			}
			K(p, j) = 0.5*K(p, j)/(h*h);
		}
	}
}