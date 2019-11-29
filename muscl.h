#ifndef MUSCL_H
#define MUSCL_H

#include <vector>

void muscl(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& ur);

#endif
