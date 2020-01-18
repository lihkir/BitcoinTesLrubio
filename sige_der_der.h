#ifndef SIGE_DER_DER
#define SIGE_DER_DER

#include "test_cases.h"

inline double sige_der_der(double phi)
{
    struct test_cases* pt_test = get_tests();
    return phi <= pt_test->phic ? 0 : pt_test->sigeddc*pow(phi/pt_test->phic, pt_test->kexp - 2);
}

#endif