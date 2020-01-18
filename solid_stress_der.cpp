#include "solid_stress_der.h"

double solid_stress_der(double phi)
{
    if (global::diff_idx == 1)
        return sige_der(phi);
    else if (global::diff_idx == 2)
        return sige_reg_der(phi);
    else
        throw std::invalid_argument("\nUndefined diff_idx!!\n"); 
}