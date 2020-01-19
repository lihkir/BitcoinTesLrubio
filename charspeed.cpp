#include <math.h> 
#include "charspeed.h"
#include "minmax_charspeed.h"
#include "test_cases.h"

double charspeed(Matrix<double>& u)
{
	struct test_cases* pt_test = get_tests();
	int M = u.cols() - 2*pt_test->gc;

	Matrix<double> *ptS = new Matrix<double>(2, M + 2*pt_test->gc); Matrix<double> &S = *ptS;
	minmax_charspeed(u, S);

	double cs = 0; double csi;
	for (int i = pt_test->gc + 1; i <= M + pt_test->gc; i++)
	{
		csi = max(abs(S(1, i)), abs(S(2, i)));
		cs = max(cs, csi);
	}

	delete ptS;

	return cs;
}