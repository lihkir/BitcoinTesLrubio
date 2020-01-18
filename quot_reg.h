#ifndef QUOT_REG
#define QUOT_REG

inline double quot_reg(double phi)
{
    struct test_cases* pt_test = get_tests();
    return pt_test->epsilon/pow(phi - pt_test->phic, 3);
}

#endif