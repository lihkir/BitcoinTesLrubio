#include <math.h>
#include <vector>
#include "diffusion_tensor.h"
#include "diffus.h"
#include "test_cases.h"
#include "utilities.h"

void diffus(std::vector<std::vector<double>>& uh, double h, std::vector<std::vector<double>> &K)
{
	struct test_cases* pt_test = get_tests();

	int m = uh.size();
	int n = uh[0].size() - 2 * pt_test->gc;

	std::vector<std::vector<double>> u = sub_matrix(uh, pt_test->gc);
	std::map<int, std::vector<std::vector<double>>> B = diffusion_tensor(u);

	for (int p = 0; p < m; p++)
	{
		for (int j = 0; j < n; j++)
		{
			K[p][j] = 0;
			if (j > 0)
			{
				for (int i = 0; i < m; i++)
				{
					K[p][j] = K[p][j] - (B[j][p][i] + B[j - 1][p][i]) * (u[i][j] - u[i][j - 1]);
				}
			}
			if (j < n - 1) {
				for (int i = 0; i < m; i++)
				{
					K[p][j] = K[p][j] + (B[j + 1][p][i] + B[j][p][i]) * (u[i][j + 1] - u[i][j]);
				}
			}
			K[p][j] = 0.5 / (h * h) * K[p][j];
		}
	}
}