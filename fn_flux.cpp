#include "fn_flux.h"

void fn_flux(Matrix<double> &phi, Matrix<double> &f)
{
	struct test_cases* pt_test = get_tests();

	int N = phi.rows();	
	Matrix<double> &delta = *pt_test->delta;
	
	double vrho = sed_hsf(sum1D(phi));
	double p2 = dot_product(phi, delta);

	for (int j = 1; j <= N; j++) f(j) = phi(j)*vrho*(delta(j) - p2);

	delete pt_test->delta;
	delete pt_test;
}