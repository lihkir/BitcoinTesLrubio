#ifndef TEST_CASES_H 
#define TEST_CASES_H

#include "containers.h"

namespace global 
{ 
	extern int idx_test; 
	extern int idx_q;
}

struct test_cases* get_tests();

struct test_cases
{
	Matrix<double> *delta;
	
	double phimax;
	double g;
	double mu_f;
	double mu_g;
	double rho_s;
	double vinf;
	double L;
	double beta;
	double nexp;
	double phi0;
	double phic;
	double sigma0;
	double kexp;
	double sigedc;
	double sigeddc;
	double epsilon;

	int convec_type;
	int M_rows;
	int gc;
	int gc_id;
};

#endif