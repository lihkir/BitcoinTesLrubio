#include <stdio.h>
#include <math.h>
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "jacobiana_dec.h"
#include "utilities.h"

std::vector<std::vector<double>> jacobiana_dec(std::vector<double>& phi)
{
	struct test_cases* pt_test = get_tests();
	int N = phi.size();

	std::vector<std::vector<double>> Jdec(N, std::vector<double>(N));

	double rho = vector_sum(phi);
	double wrho = sed_hsf(rho);
	double w1rho = sed_hsf_der(rho);

	// J=diag(SUBMATRIX(Jdec, 1))+SUBMATRIX(Jdec, 2)*SUBMATRIX(Jdec, 3)"
	for (int i = 0; i < N; i++) 
	{
		Jdec[i][0] = wrho * pt_test->delta[i];
		Jdec[i][1] = phi[i] * w1rho * pt_test->delta[i];
	}
	return Jdec;
}