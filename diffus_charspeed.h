#ifndef DIFFUS_CHARSPEED_H
#define DIFFUS_CHARSPEED_H

#include <vector>

namespace global { extern double D0; }

double diffus_charspeed(std::vector<std::vector<double>> &u);

#endif
