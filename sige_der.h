#ifndef SIGE_DER
#define SIGE_DER

#include "test_cases.h"

inline double sige_der(double phi)
{
    struct test_cases* pt_test = get_tests();
    return phi <= pt_test->phic ? 0 : pt_test->sigedc*pow(phi/pt_test->phic, pt_test->kexp - 1);
}

#endif