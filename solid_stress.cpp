#include "solid_stress.h"

double solid_stress(double phi)
{
    if (global::diff_idx == 1)
        return sige(phi);
    else if (global::diff_idx == 2)
        return sige_reg(phi);
    else
        throw std::invalid_argument("\nUndefined diff_idx!!\n"); 
}