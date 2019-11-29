#include <math.h>
#include "minmod.h"
#include "muscl.h"
#include "test_cases.h"
#include "minmod.h"

void muscl(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& ur)
{

	struct test_cases* pt_test = get_tests();
	int M = u[0].size() - 2 * pt_test->gc; int N = u.size();

	double du;
	for (int i = pt_test->gc - 1; i < M + pt_test->gc + 1; i++)
	{
		for (int j = 0; j < N; j++)
		{
			du = minmod(u[j][i] - u[j][i-1], u[j][i+1] - u[j][i]);
			ur[j][i-1] = u[j][i] - 0.5 * du;
			ul[j][i] = u[j][i] + 0.5 * du;
		}
	}
}