#ifndef SIGE_REG_DER_H
#define SIGE_REG_DER_H

#include "sige.h"
#include "sige_der.h"
#include "exp_reg.h"
#include "quot_reg.h"

inline double sige_reg_der(double phi) { return exp_reg(phi)*(sige_der(phi) + 2*sige(phi)*quot_reg(phi)); }

#endif
