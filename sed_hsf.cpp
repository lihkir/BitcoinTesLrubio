#include "sed_hsf.h"

double sed_hsf(double phi)
{
	struct test_cases* pt_test = get_tests();
	if (phi < 0 || phi > pt_test->phimax)
    	return 0;
	else if (phi > pt_test->phi0)
	{
    	double vphi0 = pow(1 - pt_test->phi0, pt_test->nexp - 2);
    	double v1phi0 = -(pt_test->nexp - 2)*pow(1 - pt_test->phi0, pt_test->nexp - 3);
    	return pt_test->vinf*(1 - phi)*(vphi0 + v1phi0*(phi - pt_test->phi0)); 
	}
	else 
		return pt_test->vinf*pow(1 - phi, pt_test->nexp - 1);
}