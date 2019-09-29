#ifndef SED_HSF_DER_H
#define SED_HSF_DER_H

#include <cmath>
#include "test_cases.h"

inline double sed_hsf_der(double phi)
{
	struct test_cases* pt_test = get_tests();
	return phi >= 0 && phi <= 1 ? -pt_test->nexp * pt_test->vinf * pow(1 - phi, pt_test->nexp - 1) : 0;
}

#endif
