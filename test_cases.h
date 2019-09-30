#ifndef TEST_CASES_H 
#define TEST_CASES_H

#include <vector>

namespace global { extern int idx_test; }
struct test_cases* get_tests();

struct test_cases
{
	std::vector<double> delta;
	int M_rows = NULL;
	int gc = NULL;
	double phimax = NULL;
	int int_form = NULL;
	int gc_id = NULL;
	double g = NULL;
	double mu_c = NULL;
	double rho_d = NULL;
	double rho_c = NULL;
	double vinf = NULL;
	double L = NULL;
	double nexp = NULL;
};

#endif