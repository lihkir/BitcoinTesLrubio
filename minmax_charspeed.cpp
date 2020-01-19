#include "minmax_charspeed.h"
#include "update_utilities.h"

void minmax_charspeed(Matrix<double> &ua, Matrix<double> &Sa)
{
	struct test_cases* pt_test = get_tests();
	
	int N2 = ua.rows();
	int M2 = ua.cols() - 2 * pt_test->gc;

	Matrix<double> *ptul = new Matrix<double>(N2, 1); Matrix<double> &ul = *ptul;
	Matrix<double> *ptur = new Matrix<double>(N2, 1); Matrix<double> &ur = *ptur;
	Matrix<double> *ptuh = new Matrix<double>(3, 1);  Matrix<double> &uh = *ptuh;
	Matrix<double> *ptSh = new Matrix<double>(2, M2 + 2 * pt_test->gc); Matrix<double> &Sh = *ptSh;
	
	for (int i = pt_test->gc + 1; i <= M2 + pt_test->gc; i++)
	{
		get_col(ul, ua, i);
		get_col(ur, ua, i + 1);
		lminmax_charspeed(ul, ur, uh);
		for (int j = 1; j <= 2; j++) Sa(j, i) = uh(j);
	}
	delete ptSh;
	delete ptuh;
	delete ptul;
	delete ptur;
}