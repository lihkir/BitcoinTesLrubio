#include <math.h>
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "jacobiana.h"
#include "utilities.h"

std::vector<std::vector<double>> jacobiana(std::vector<double>& phi)
{
	struct test_cases* pt_test = get_tests();
	int N = phi.size();

	std::vector<std::vector<double>> jac(N, std::vector<double>(N));

	double rho = VectorSum(phi);
	double u = sed_hsf(rho);
	double v = sed_hsf_der(rho);

	for (int i = 0; i < N; i++) 
		for (int j = 0; j < N; j++) 
			jac[i][j] = phi[i] * v * pt_test->delta[i];

	for (int i = 0; i < N; i++) jac[i][i] = jac[i][i] + u * pt_test->delta[i];
	
	return jac;
}