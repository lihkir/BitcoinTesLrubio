#ifndef SIGE_H
#define SIGE_H

#include "test_cases.h"

inline double sige(double phi)
{
    struct test_cases* pt_test = get_tests();
    return phi <= pt_test->phic ? 0 : pt_test->sigma0*(pow(phi/pt_test->phic, pt_test->kexp) - 1);
}

#endif
