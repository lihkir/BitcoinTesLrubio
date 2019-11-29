#include <stdio.h>
#include "update_utilities.h"

void UpdateCol(std::vector<std::vector<double>> &m, std::vector<std::vector<double>> &c, int col)
{
    for (unsigned int i = 0; i < m.size(); i++) m[i][col] = c[i][0];
}