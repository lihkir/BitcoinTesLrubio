#ifndef _DIFFUS_CHARSPEED_H_
#define _DIFFUS_CHARSPEED_H_

#include <vector>

namespace global { extern double D0; }

double diffus_charspeed(std::vector<std::vector<double>> &u);

#endif
