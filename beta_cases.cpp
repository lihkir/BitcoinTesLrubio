#include "beta_cases.h"

std::vector<double> beta_cases(int test_id, int& M_rows, std::vector<double>* pt_beta)
{
	if (test_id == 1)
	{
		pt_beta = new std::vector<double>();
		pt_beta->push_back(60.0);
		pt_beta->push_back(67.5);
		pt_beta->push_back(75.0);
		pt_beta->push_back(82.5);
		pt_beta->push_back(90.0);
		pt_beta->push_back(97.5);
		pt_beta->push_back(105.0);
		pt_beta->push_back(112.5);
		// pt_beta->push_back(120.0);
		M_rows = pt_beta->size();
	}
	std::vector<double>& beta = *pt_beta;
	return beta;
}