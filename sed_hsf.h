#ifndef SED_HSF
#define SED_HSF

#include <cmath>
#include "test_cases.h"

inline double sed_hsf(double phi)
{
	struct test_cases* pt_test = get_tests();
	return exp(-0.5 * pow(phi / pt_test->phimax, 2));
}

#endif