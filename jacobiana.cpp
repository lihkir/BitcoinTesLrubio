#include "jacobiana.h"
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "test_cases.h"

void jacobiana(Matrix<double> &phi, Matrix<double> &jac)
{
	struct test_cases* pt_test = get_tests();

	Matrix<double> &delta = *pt_test->delta;
	int N = size(phi, 1);

	double rho = sum1D(phi);
	double p2 = dot_product(phi, delta);

	double u = sed_hsf(rho);
	double v = sed_hsf_der(rho); 

	for (int i = 1; i <= N; i++)
		for (int j = 1; j <= N; j++)
			jac(i, j) = phi(i)*(v*(delta(i) - p2) - u*delta(j));

	for (int i = 1; i <= N; i++)	
    	jac(i, i) = jac(i, i) + u*(delta(i) - p2);

	delete pt_test->delta;
	delete pt_test;
}