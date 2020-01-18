#ifndef SIGE_REG_DER_DER_H
#define SIGE_REG_DER_DER_H

#include "sige.h"
#include "sige_der.h"
#include "sige_der_der.h"
#include "exp_reg.h"
#include "quot_reg.h"

inline double sige_reg_der_der(double phi) 
{ 
    struct test_cases* pt_test = get_tests();
    return exp_reg(phi)*(sige_der_der(phi) + 4*sige_der(phi)*quot_reg(phi) + 2*pt_test->epsilon*sige(phi)*(2*pt_test->epsilon - 3*pow(phi - pt_test->phic, 2))); 
}

#endif
