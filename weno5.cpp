#include <stdio.h>
#include <math.h>
#include <stdexcept>
#include "weno5.h"

double weno5(std::vector<double>& f, double eps)
{
	double I0 = f[0] * (4.0 * f[0] - 19.0 * f[1] + 11.0 * f[2]) + f[1] * (25.0 * f[1] - 31.0 * f[2]) + f[2] * (10.0 * f[2]);
	I0 = I0 / 3.0;

	double p0 = 1.0 / 3.0 * f[0] - 7.0 / 6.0 * f[1] + 11.0 / 6.0 * f[2];

	double I1 = f[1] * (4.0 * f[1] - 13.0 * f[2] + 5.0 * f[3]) + f[2] * (13.0 * f[2] - 13.0 * f[3]) + f[3] * (4.0 * f[3]);
	I1 = I1 / 3.0;

	double p1 = -1.0 / 6.0 * f[1] + 5.0 / 6.0 * f[2] + 1.0 / 3.0 * f[3];

	double I2 = f[2] * (10.0 * f[2] - 31.0 * f[3] + 11.0 * f[4]) + f[3] * (25.0 * f[3] - 19.0 * f[4]) + f[4] * (4.0 * f[4]);
	I2 = I2 / 3.0;

	double p2 = 1.0 / 3.0 * f[2] + 5.0 / 6.0 * f[3] - 1.0 / 6.0 * f[4];

	double a0 = (I0 + eps) * (I0 + eps);
	a0 = 1.0 / 10.0 / a0;
	double a1 = (I1 + eps) * (I1 + eps);
	a1 = 3.0 / 5.0 / a1;
	double a2 = (I2 + eps) * (I2 + eps);
	a2 = 3.0 / 10.0 / a2;

	double sum_a = a0 + a1 + a2;
	double w0 = a0 / sum_a;
	double w1 = a1 / sum_a;
	double w2 = a2 / sum_a;
	double w = w0 * p0 + w1 * p1 + w2 * p2;

	//printf("WENO: % .15e % .15e % .15e % .15e % .15e % .15e\n", f[0], f[1], f[2], f[3], f[4], w);

	if (isnan(w)) throw std::invalid_argument("\nisnan(w)!!\n");

	return w;
}