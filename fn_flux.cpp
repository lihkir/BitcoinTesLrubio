#include "fn_flux.h"
#include "sed_hsf.h"
#include "utilities.h"
#include "test_cases.h"

std::vector<double> fn_flux(std::vector<double> phi)
{
	struct test_cases* pt_test = get_tests();
	std::vector<double> f(phi.size());

	double phit = vector_sum(phi);
	double vphi = sed_hsf(phit);
	double dphi = dot_product(pt_test->delta, phi);

	for (unsigned int j = 0; j < phi.size(); j++)
		f[j] = pt_test->mudrho*phi[j]*vphi*(1 - phit)*(pt_test->delta[j] - dphi);
	
	return f;
}