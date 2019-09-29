#include <stdio.h>
#include <math.h>
#include "lminmax_charspeed.h"
#include "minmax_charspeed.h"
#include "test_cases.h"

std::vector<std::vector<double>> minmax_charspeed(std::vector<std::vector<double>>& ua)
{
	struct test_cases* pt_test = get_tests();
	int N2 = ua.size();
	int M2 = ua[0].size() - 2 * pt_test->gc;

	std::vector<double> ul(N2), ur(N2), uh(N2);
	std::vector<std::vector<double>> Sh(2, std::vector<double>(M2 + 2 * pt_test->gc));

	for (int i = pt_test->gc; i < M2 + pt_test->gc; i++)
	{
		for (int j = 0; j < N2; j++)
		{
			ul[j] = ua[j][i];
			ur[j] = ua[j][i + 1];
		}
		std::vector<double> uh = lminmax_charspeed(ul, ur);
		for (int j = 0; j < 2; j++) Sh[j][i] = uh[j];
	}
	return Sh;
}