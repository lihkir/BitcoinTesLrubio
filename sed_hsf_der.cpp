#include "sed_hsf_der.h"

double sed_hsf_der(double phi)
{
	struct test_cases* pt_test = get_tests();

    if (phi < 0 || phi > pt_test->phimax)
        return 0;
    else if (phi > pt_test->phi0)
    {
        double vphi0 = pow(1 - pt_test->phi0, pt_test->nexp - 2);
        double v1phi0 = -(pt_test->nexp - 2)*pow(1 - pt_test->phi0, pt_test->nexp - 3);
        return pt_test->vinf*(-vphi0 + v1phi0*(1 + pt_test->phi0) - 2*phi*v1phi0);
    }
    else
        return -(pt_test->nexp - 1)*pt_test->vinf*pow(1 - phi, pt_test->nexp - 2);
}