#include <iostream>
#include <cstdio>
#include <cmath>
#include "test_cases.h"

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
		pt_test->gc_id = 1;		
		pt_test->L = 0.02;
		pt_test->mu_f = 1e-3;
		pt_test->rho_s = 1800;
		pt_test->g = 9.81;
		pt_test->phimax = 0.66;	
		pt_test->nexp = 4.65;
		pt_test->sigma0 = 180;
		pt_test->k = 2;
		pt_test->phic = 0.2;
	  	pt_test->phi0 = ((pt_test->nexp-2)*pt_test->phimax-1)/(pt_test->nexp-3);
		pt_test->vphis = pow((1 - pt_test->phi0), (pt_test->nexp - 2));
        pt_test->vphis_der = -(pt_test->nexp - 2)*pow((1 - pt_test->phi0), (pt_test->nexp - 3));
		pt_test->mudrho = pt_test->mu_f*pt_test->rho_s;
		pt_test->muog = pt_test->mu_f/pt_test->g;
		
		pt_test->delta.push_back(1);
		pt_test->delta.push_back(0.5);
		pt_test->delta.push_back(0.25);
		pt_test->M_rows = pt_test->delta.size();
	}
	else if (global::idx_test == 2 ) 
	{
		pt_test->gc_id = 1;		
		pt_test->L = 0.02;
		pt_test->mu_f = 1e-3;
		pt_test->rho_s = 1800;
		pt_test->g = 9.81;
		pt_test->phimax = 0.66;	
		pt_test->nexp = 4.65;
		pt_test->sigma0 = 180;
		pt_test->k = 2;
		pt_test->phic = 0.2;
	  	pt_test->phi0 = ((pt_test->nexp-2)*pt_test->phimax-1)/(pt_test->nexp-3);
		pt_test->vphis = pow((1 - pt_test->phi0), (pt_test->nexp - 2));
        pt_test->vphis_der = -(pt_test->nexp - 2)*pow((1 - pt_test->phi0), (pt_test->nexp - 3));

		pt_test->delta.push_back(1);
		pt_test->delta.push_back(0.5);
		pt_test->delta.push_back(0.25);
		pt_test->M_rows = pt_test->delta.size();   
	}
	else if (global::idx_test == 3) 
	{
		pt_test->gc_id = 1;		
		pt_test->L = 0.02;
		pt_test->mu_f = 1e-3;
		pt_test->rho_s = 1800;
		pt_test->g = 9.81;
		pt_test->phimax = 0.66;	
		pt_test->nexp = 4.65;
		pt_test->sigma0 = 180;
		pt_test->k = 2;
		pt_test->phic = 0.2;
	  	pt_test->phi0 = ((pt_test->nexp-2)*pt_test->phimax-1)/(pt_test->nexp-3);
		pt_test->vphis = pow((1 - pt_test->phi0), (pt_test->nexp - 2));
        pt_test->vphis_der = -(pt_test->nexp - 2)*pow((1 - pt_test->phi0), (pt_test->nexp - 3));

		pt_test->delta.push_back(1);
		pt_test->delta.push_back(0.5);
		pt_test->delta.push_back(0.25);
		pt_test->M_rows = pt_test->delta.size();
	}
	else if (global::idx_test == 4) 
	{
		pt_test->gc_id = 1;		
		pt_test->L = 0.02;
		pt_test->mu_f = 1e-3;
		pt_test->rho_s = 1800;
		pt_test->g = 9.81;
		pt_test->phimax = 0.66;	
		pt_test->nexp = 4.65;
		pt_test->sigma0 = 180;
		pt_test->k = 2;
		pt_test->phic = 0.2;
	  	pt_test->phi0 = ((pt_test->nexp-2)*pt_test->phimax-1)/(pt_test->nexp-3);
		pt_test->vphis = pow((1 - pt_test->phi0), (pt_test->nexp - 2));
        pt_test->vphis_der = -(pt_test->nexp - 2)*pow((1 - pt_test->phi0), (pt_test->nexp - 3));

		pt_test->delta.push_back(1);
		pt_test->delta.push_back(0.5);
		pt_test->delta.push_back(0.25);
		pt_test->M_rows = pt_test->delta.size();
	}
	else 
	{
		std::cerr<<"Test "<<global::idx_test<<" no definido"<<std::endl;
		exit(1);
	}
	return pt_test;
}