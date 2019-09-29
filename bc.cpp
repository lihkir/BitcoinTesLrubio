#include "bc.h"
#include "test_cases.h"

void bc(std::vector<std::vector<double>>& u)
{
	struct test_cases* pt_test = get_tests();
	int N = u[0].size() - 2 * pt_test->gc;

	for (unsigned int i = 0; i < u.size(); i++)
	{
		for (auto j = 0; j != pt_test->gc; j++)
		{
			if (pt_test->gc_id == 1)
			{
				u[i][j] = u[i][N + j];
				u[i][pt_test->gc + N + j] = u[i][pt_test->gc + j];
			}
			else
			{
				u[i][j] = 1e10;
				u[i][pt_test->gc + N + j] = 1e10;
			}
		}
	}
}