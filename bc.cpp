#include "bc.h"

void bc(Matrix<double>& u)
{
	struct test_cases* pt_test = get_tests();

	int m = u.rows();
	int n = u.cols() - 2*pt_test->gc;

	for (int j = 1; j <= pt_test->gc; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			u(i, j) = u(i, 2*pt_test->gc + 1 - j);
			u(i, pt_test->gc + n + j) = u(i, pt_test->gc + 1 + n - j);			
		}		
	}
}