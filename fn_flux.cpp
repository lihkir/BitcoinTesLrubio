#include "fn_flux.h"
#include "sed_hsf.h"
#include "utilities.h"
#include "test_cases.h"

std::vector<double> fn_flux(std::vector<double> phi)
{
	struct test_cases* pt_test = get_tests();
	std::vector<double> f(phi.size());
	double vrho = sed_hsf(VectorSum(phi));
	for (unsigned int j = 0; j < phi.size(); j++) f[j] = phi[j] * vrho * pt_test->delta[j];
	return f;
}