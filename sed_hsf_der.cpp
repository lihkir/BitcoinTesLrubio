#include <cmath>
#include "stdio.h"
#include "sed_hsf_der.h"
#include "test_cases.h"

double sed_hsf_der(double phi)
{
	struct test_cases* pt_test = get_tests();

	if (phi < 0 || phi > pt_test->phimax) 
        return 0;
    else if (phi > pt_test->phi0)
        return pt_test->vphis_der;
    else 
        return -(pt_test->nexp-2)*pow((1 - phi), (pt_test->nexp - 3));
}
