#include "test_cases.h"
#include <iostream>
#include <cstdio>

struct test_cases* get_tests()
{
	struct test_cases* pt_test = new test_cases;

	if (global::idx_q == -1) {
		pt_test->convec_type = 1; pt_test->gc = 3;
	} else {
		pt_test->convec_type = 0; pt_test->gc = 2;
	}

	if (global::idx_test == 1)
	{
		pt_test->M_rows = 3;
		pt_test->gc_id = 1;
		pt_test->g = 9.81;
		pt_test->mu_f = 1e-3; 
		pt_test->rho_s = 1800;
		pt_test->beta = 1e-15;
		double d1 = 1.19e-5;
		pt_test->vinf = pt_test->g*d1*d1*pt_test->rho_s/(18*pt_test->mu_f);
		pt_test->mu_g = pt_test->vinf/pt_test->rho_s;
		pt_test->L = 1;
		pt_test->phimax = 0.66;
		pt_test->nexp = 4.7;
		pt_test->phi0 = ((pt_test->nexp - 2)*pt_test->phimax-1)/(pt_test->nexp - 3);
		pt_test->phic = 0.2;
		pt_test->sigma0 = 180;
		pt_test->kexp = 2;
		pt_test->sigedc = pt_test->sigma0*pt_test->kexp/pt_test->phic;
		pt_test->sigeddc = pt_test->sigma0*pt_test->kexp*(pt_test->kexp - 1)/pow(pt_test->phic, 2);
		pt_test->epsilon = 1e-5;
		
		pt_test->delta =  new Matrix<double>(pt_test->M_rows, 1); Matrix<double> &delta = *pt_test->delta;
		delta(1) = 1;
		delta(2) = 0.5;
		delta(3) = 0.25;
	}
	else 
	{
		std::cerr<<"Test "<<global::idx_test<<" no definido"<<std::endl;
		exit(1);
	}
	return pt_test;
}