#include <stdio.h>
#include <math.h>
#include "minmax_charspeed.h"
#include "charspeed.h"
#include "test_cases.h"
#include "utilities.h"

double charspeed(std::vector<std::vector<double>>& u)
{
	struct test_cases* pt_test = get_tests();
	auto gc2 = 2 * pt_test->gc; auto M = u[0].size() - gc2;

	double csm = 0;
	std::vector<std::vector<double>> Sminmax = minmax_charspeed(u);
	for (auto i = pt_test->gc; i < M + pt_test->gc - 1; i++)
		csm = max(csm, max(fabs(Sminmax[0][i]), fabs(Sminmax[1][i])));

	return csm;
}

