#ifndef SOLID_STRESS_DER_DER_H
#define SOLID_STRESS_DER_DER_H

#include "sige_der_der.h"
#include "sige_reg_der_der.h"

namespace global { extern int diff_idx; }

double solid_stress_der_der(double phi);

#endif