#ifndef SIGE_REG_H
#define SIGE_REG_H

#include "sige.h"
#include "exp_reg.h"

inline double sige_reg(double phi) { return sige(phi)*exp_reg(phi); }

#endif