#include "fn_flux.h"
#include "sed_hsf.h"
#include "utilities.h"
#include <stdio.h>

std::vector<double> fn_flux(std::vector<double> phi, std::vector<double> beta)
{
	std::vector<double> f(phi.size());
	double vrho = sed_hsf(VectorSum(phi));
	for (unsigned int j = 0; j < phi.size(); j++) f[j] = phi[j] * vrho * beta[j];
	return f;
}
