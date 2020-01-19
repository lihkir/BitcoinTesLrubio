#include <math.h>
#include "minmod.h"
#include "muscl.h"
#include "test_cases.h"
#include "minmod.h"

void muscl(Matrix<double>& u, Matrix<double>& ul, Matrix<double>& ur)
{
	struct test_cases* pt_test = get_tests();
	
	int M = size(u, 2) - 2*pt_test->gc;
	int N = size(u, 1);

	double du;
	for (int i = pt_test->gc; i <= M + pt_test->gc + 1; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			du = minmod(u(j, i)-u(j, i-1), u(j, i+1)-u(j, i));
    		ur(j, i-1)=u(j, i)-0.5*du;
    		ul(j, i)=u(j, i)+0.5*du;
		}
	}
}