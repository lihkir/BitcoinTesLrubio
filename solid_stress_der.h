#ifndef SOLID_STRESS_DER_H
#define SOLID_STRESS_DER_H

#include "test_cases.h"
#include <cmath>

inline double solid_stress_der(double phi)
{
    struct test_cases *pt_test = get_tests();
    return phi <= pt_test->phic ? 0 : (pt_test->sigma0*pt_test->k/pt_test->phic)*pow(phi/pt_test->phic, pt_test->k-1);
}

#endif