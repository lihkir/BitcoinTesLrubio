#include <cmath>
#include "stdio.h"
#include "sed_hsf.h"
#include "test_cases.h"

double sed_hsf(double phi)
{
	struct test_cases* pt_test = get_tests();

    if (phi < 0 || phi > pt_test->phimax)
        return 0;
    else if (phi > pt_test->phi0)
        return pt_test->vphis + pt_test->vphis_der*(phi - pt_test->phi0);
    else
        return pow((1 - phi), (pt_test->nexp - 2));
}
