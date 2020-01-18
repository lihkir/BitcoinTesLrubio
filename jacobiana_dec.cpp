#include <stdio.h>
#include <math.h>
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "jacobiana_dec.h"

void jacobiana_dec(Matrix<double> &phi, Matrix<double> &Jdec)
{
	struct test_cases* pt_test = get_tests();
	
	Matrix<double> &delta = *pt_test->delta;
	int N = size(phi, 1);
	
	double rho = sum1D(phi);
	double p2 = dot_product(phi, delta);

	double wrho = sed_hsf(rho);
	double w1rho = sed_hsf_der(rho);

	/** J=diag(Jdec(:,1))+Jdec(:, 2:3)*Jdec(:, 4:5) **/
	for (int i = 1; i <= N; i++)
	{
		Jdec(i, 1) = wrho*(delta(i)-p2);
		Jdec(i, 2) = phi(i)*w1rho*(delta(i)-p2);
  		Jdec(i, 3) = phi(i)*wrho;
  		Jdec(i, 4) = 1;
  		Jdec(i, 5) = -delta(i);
	}
}