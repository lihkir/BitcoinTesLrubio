#ifndef EXP_REG
#define EXP_REG

inline double exp_reg(double phi)
{
    struct test_cases* pt_test = get_tests();
    return exp(-pt_test->epsilon/pow(phi - pt_test->phic, 2));
}

#endif