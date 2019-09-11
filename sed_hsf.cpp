#include<cmath>
#include "sed_hsf.h"
#include "globals.h"

double sed_hsf(double phi)
{
  global_data* pt_data;
  return exp(-0.5 * pow(phi / pt_data->phimax, 2));
}
