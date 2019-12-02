#ifndef TEST_CASES_H 
#define TEST_CASES_H

#include <vector>

namespace global 
{ 
	extern int idx_test; 
	extern int idx_q;
}

struct test_cases* get_tests();

struct test_cases
{
	std::vector<double> delta;
	int M_rows;
	int gc;
	double phimax;
	int gc_id;
	double g;
	double L;
	double nexp;
	int convec_type;
	double sigma0;
	double k;
	double phic;
	double mu_f;
	double rho_s;
	double phi0;
	double vphis;
	double vphis_der;
	double mudrho;
	double muog;
};

#endif