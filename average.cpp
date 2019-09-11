#include "average.h"
#include <stdio.h>

std::vector<double> average(std::vector<double> ul, std::vector<double> ur, double s)
{
	std::vector<double> ua(ul.size());
	for (unsigned int i = 0; i < ua.size(); i++)
		ua[i] = ul[i] + s*(ur[i]-ul[i]);

	return ua;
}
