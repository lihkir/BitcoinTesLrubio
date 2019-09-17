#include <stdio.h>
#include <math.h>
#include <vector>
#include "lminmax_charspeed.h"
#include "minmax_charspeed.h"
#include "test_cases.h"
#include "utilities.h"

void minmax_charspeed(std::vector<std::vector<double>> &ua, std::vector<std::vector<double>> &Sa)
{
    struct test_cases *pt_test = get_tests();
    int N2 = ua.size();
    int M2 = ua[0].size() - 2*pt_test->gc;

    std::vector<double> ul(N2), ur(N2), uh(N2);

    for (int i = pt_test->gc; i < M2 + pt_test->gc; i++) 
    {
        printf("\n\ni=%d\t M2 = %d\t N2 = %d\n\n", i, M2, N2);
        for (int j = 0; j < N2; j++) 
        {
            ul[j] = ua[j][i];
            ur[j] = ua[j][i+1];
        }
        PrintingContainer(ul);
        PrintingContainer(ur);
        PrintingContainer(ua);
        std::vector<double> uh = lminmax_charspeed(ul, ur);
        PrintingContainer(uh);
        for (int j = 0; j < 2; j++) 
        {
            Sa[j][i] = uh[j];
        }
        PrintingContainer(Sa);
    }
}