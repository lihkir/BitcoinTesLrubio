#ifndef TEST_CASES_H 
#define TEST_CASES_H

#include <vector>

struct test_cases *get_tests();

struct test_cases
{
	std::vector<double> beta;
    std::vector<double> delta;
    int M_rows;
    int gc;
    double phimax;
    int int_form;
    int gc_id;
    double g;
    double mu_c;
    double rho_d;
    double rho_c;
    double vinf;
    double L;
    double nexp;
};

#endif