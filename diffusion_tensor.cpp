#include <stdio.h>
#include <math.h>
#include <assert.h> 
#include <map>
#include "test_cases.h"
#include "utilities.h"
#include "diffusion_tensor.h"
#include "sed_hsf.h"
#include "solid_stress.h"
#include "solid_stress_der.h"
#include "deltak.h"

std::map<int, std::vector<std::vector<double>>> diffusion_tensor(std::vector<std::vector<double>>& u)
{
	struct test_cases* pt_test = get_tests();

	int m = u.size();
	int n = u[0].size();

	std::map<int, std::vector<std::vector<double>>> Bsol;
	std::vector<std::vector<double>> Bele(m, std::vector<double>(m));

	for (int k = 0; k < n; k++)
	{	
		std::vector<double> phi = col(u, k); 
		double dphi = dot_product(pt_test->delta, phi);
		double phit = vector_sum(phi);
		double vphi = sed_hsf(phit);
		double vopq = vphi/phit;
		double sige = solid_stress(phit);
		double sigd = solid_stress_der(phit);

		if (phit > 1) throw std::invalid_argument("\nError: phit > 1 !!\n");
		Bsol[k] = Bele;
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				Bsol[k][i][j] = pt_test->muog*vopq*((1 - phit)*phi[i]*(pt_test->delta[i] - dphi)*sigd - (pt_test->delta[i]*deltak(i, j) - pt_test->delta[j]*phi[i] - (phi[i]/phit)*(pt_test->delta[i]-dphi))*sige);
	}
	return Bsol;
}