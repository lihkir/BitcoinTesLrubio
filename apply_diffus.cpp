#include <stdio.h>
#include <math.h>
#include <vector>
#include "diffusion_tensor.h"
#include "apply_diffus.h"
#include "test_cases.h"
#include "utilities.h"

std::vector<std::vector<double>> apply_diffus(std::vector<std::vector<double>> &ut, std::vector<std::vector<double>> &uh, double h)
{
	struct test_cases* pt_test = get_tests();
	
	auto m = uh.size();
	auto n = (int)uh[0].size() - 2 * pt_test->gc;
	std::vector<std::vector<double>> K(m, std::vector<double>(n));

	std::vector<std::vector<double>> u   = SubMatrix(uh);
	std::vector<std::vector<double>> ut0 = SubMatrix(ut);

	std::map<int, std::vector<std::vector<double>>> B = diffusion_tensor(ut0);

	for (auto p = 0; p < m; p++) 
	{
		for (auto j = 0; j < n; j++) 
		{
			K[p][j] = 0;
			if (j > 0) 
			{
				for (auto i = 0; i < m; i++) 
				{
					K[p][j] = K[p][j] - (B[j][p][i] + B[j-1][p][i]) * (u[i][j] - u[i][j-1]);
				}
			}
			if (j < n-1) {
				for (auto i = 0; i < m; i++) 
				{
					K[p][j] = K[p][j] + (B[j+1][p][i] + B[j][p][i]) * (u[i][j+1] - u[i][j]);
				}
			}
			K[p][j] = 0.5 / (h * h) * K[p][j];
		}
	}
	return K;
}