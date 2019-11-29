#include <stdio.h>
#include <math.h>
#include <assert.h> 
#include <map>
#include "diffusion_tensor.h"
#include "test_cases.h"
#include "utilities.h"

std::map<int, std::vector<std::vector<double>>> diffusion_tensor(std::vector<std::vector<double>>& u)
{
	struct test_cases* pt_test = get_tests();

	int m = u.size();
	int n = u[0].size();

	std::map<int, std::vector<std::vector<double>>> Bsol;
	std::vector<std::vector<double>> Bele(m, std::vector<double>(m));

	for (int k = 0; k < n; k++)
	{
		double phit = VectorSum(Col(u, k));
		if (phit > 1)
		{	
			printf("\nInside diffusion_tensor phit = %0.16f\n", phit);
			throw std::invalid_argument("\nError: phit > 1 !!\n");
		}

		Bsol[k] = Bele;
		for (auto i = 0; i < m; i++)
			Bsol[k][i][i] = global::D0 * pow(1 - phit, pt_test->nexp);
	}
	return Bsol;
}