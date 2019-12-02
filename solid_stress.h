#ifndef SOLID_STRESS_H
#define SOLID_STRESS_H

#include "test_cases.h"
#include <cmath>

inline double solid_stress(double phi)
{
    struct test_cases *pt_test = get_tests();
    return phi <= pt_test->phic ? 0 : pt_test->sigma0*(pow(phi/pt_test->phic, pt_test->k) - 1);
}

#endif