#include <math.h>
#include "minmod.h"
#include "utilities.h"

double minmod(double a, double b)
{
	double m;

	if (a * b < 0) 
		m = 0;
	else 
	{
		m = min(fabs(a), fabs(b));
		if (a < 0)
			m = -m;
	}
	return m;
}