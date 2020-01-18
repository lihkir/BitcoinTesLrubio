#include <stdio.h>
#include "average.h"

void average(Matrix<double> &ul, Matrix<double> &ur, double s, Matrix<double> &ua)
{
	int N = size(ul, 1);
	for (int i = 1; i <= N; i++) ua(i) = ul(i) + s * (ur(i) - ul(i));
}