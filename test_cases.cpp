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
		pt_test->phimax = 50;
		pt_test->gc_id = 1;
		pt_test->g = 9.81;
		pt_test->mu_c = 6.5e-3;
		pt_test->rho_d = 1090;
		pt_test->rho_c = 880;
		pt_test->vinf = pt_test->g * (pt_test->rho_d - pt_test->rho_c) / (18 * pt_test->mu_c);
		pt_test->L = 0.02;
		pt_test->nexp = 4.65;
		
		pt_test->delta.push_back(4.057404490000000e-08);
		pt_test->delta.push_back(1.965604000000000e-08);
		pt_test->delta.push_back(9.950262001000001e-09);
		pt_test->delta.push_back(4.759068196000001e-09);
		pt_test->delta.push_back(2.341688881000000e-09);
		pt_test->delta.push_back(1.168614225000000e-09);
		pt_test->delta.push_back(5.669161000000000e-10);
		pt_test->delta.push_back(3.722342121000000e-11);
		pt_test->M_rows = pt_test->delta.size();
	}
	else if (global::idx_test == 2 ) 
	{
		pt_test->g = 9.81;
   		pt_test->mu_c = 6.5e-3;
   		pt_test->rho_d = 1090;
   		pt_test->rho_c = 880;
   		double d2=1.0000e-08;
		pt_test->vinf = d2*pt_test->g * (pt_test->rho_d - pt_test->rho_c) / (18 * pt_test->mu_c);
   		pt_test->L = 1;
   		pt_test->nexp = 4.65;

   		pt_test->delta.push_back(1);
   		pt_test->delta.push_back(6.4000e-01);
   		pt_test->delta.push_back(3.6000e-01);
		pt_test->M_rows = pt_test->delta.size();   
	}
	else if (global::idx_test == 3) 
	{
		pt_test->g = 9.81;
		pt_test->mu_c = 6.5e-3;
		pt_test->rho_d = 1090;
   		pt_test->rho_c = 880;
   		pt_test->vinf = pt_test->g * (pt_test->rho_d - pt_test->rho_c) / (18 * pt_test->mu_c);
   		pt_test->L = 0.02;
   		pt_test->nexp = 4.65;
		
		pt_test->delta.push_back(1.745725584039578e-07);
   		pt_test->delta.push_back(8.502483987726233e-08);
   		pt_test->delta.push_back(4.114961917286227e-08);
   		pt_test->delta.push_back(2.055894601643213e-08);
   		pt_test->delta.push_back(1.002369639342369e-08);
   		pt_test->delta.push_back(4.709986484700493e-09);
   		pt_test->delta.push_back(2.328921953264102e-09);
   		pt_test->delta.push_back(1.148277502313191e-09);
		pt_test->M_rows = pt_test->delta.size();
	}
	else if (global::idx_test == 4) 
	{
		pt_test->g = 9.81;
   		pt_test->mu_c = 6.5e-3;
   		pt_test->rho_d = 1090;
   		pt_test->rho_c = 880;
   		pt_test->vinf = pt_test->g * (pt_test->rho_d - pt_test->rho_c) / (18 * pt_test->mu_c);
   		pt_test->L = 0.02;
   		pt_test->nexp = 4.65;
		
		pt_test->delta.push_back(5e-5);
		pt_test->M_rows = pt_test->delta.size();
	}
	else 
	{
		std::cerr<<"Test "<<global::idx_test<<" no definido"<<std::endl;
		exit(1);
	}
	return pt_test;
}