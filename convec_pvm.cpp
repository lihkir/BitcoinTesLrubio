#include "bc.h"
#include "minmax_charspeed.h"
#include "muscl.h"
#include "numflux_hll.h"
#include "numflux_pvm.h"
#include "convec_pvm.h"
#include "fhroe.h"
#include "update_utilities.h"

void convec_pvm(Matrix<double>& v, double h, Matrix<double>& K)
{
	struct test_cases* pt_test = get_tests();

	int M = size(v, 2) - 2*pt_test->gc;
	int N = size(v, 1);

	int n_coef;
	if (global::idx_q == 2)
		n_coef = 1;
	else if (global::idx_q == 3)
		n_coef = 2;
	else if (global::idx_q == 4)
		n_coef = 3;
	else if (global::idx_q == 5)
		n_coef = 3;
	else if (global::idx_q == 6)
		n_coef = 5;
	else
		n_coef = 5;
	
	Matrix<double> S(2, M + 2*pt_test->gc);
	Matrix<double> ul(N, M + 2 * pt_test->gc);
	Matrix<double> ur(N, M + 2 * pt_test->gc);
	Matrix<double> fh(N, M + 2 * pt_test->gc);
	Matrix<double> fh1(N, 1);

	Matrix<double> *ptA;
	int N5;
	if (N == 1)
		ptA = new Matrix<double>(N, N);
	else
	{
		N5 = max(N, 5);
		ptA = new Matrix<double>(N, N5);
	}
	Matrix<double> &A = *ptA;

	Matrix<double> Al(N, N); 
	Matrix<double> fl(N, 1);
	Matrix<double> fr(N, 1);
	Matrix<double> ua(N, 1);
	Matrix<double> dif(N, 1);
	Matrix<double> uli(N, 1);
	Matrix<double> uri(N, 1);
	Matrix<double> Qurmul(N, 1);
	Matrix<double> Qc(n_coef, 1);
	Matrix<double> Si(2, 1);
	Matrix<double> Sl(2, 1);
	Matrix<double> Sr(2, 1);

	bc(v);
	minmax_charspeed(v, S);

	double cs=0;
	double csi;
	for (int i = pt_test->gc + 1; i <= M + pt_test->gc; i++)
	{
		csi = max(abs(S(1, i)), abs(S(2, i)));
  		cs = max(cs, csi);
	}

	muscl(v, ul ,ur);

	for (int i = pt_test->gc + 1; i <= M + pt_test->gc - 1; i++)
	{
		get_col(uli, ul, i);
		get_col(uri, ur, i);
		get_col(Sl, S, i);
		get_col(Sr, S, i + 1);

		if (global::idx_q > 0) {
			numflux_pvm(uli, uri, Sl, Sr, global::idx_q, fh1, A, Al, fl, fr, ua, dif, Qurmul, Qc, Si);
		} else if (global::idx_q == -2) {
			fhroe(uli, uri, fl, fr, fh);
		} else {
			numflux_hll(uli, uri, Sl, Sr, fh1, fl, fr, Si);
		}
		for (int j = 1; j <= N; j++) fh(j, i) = fh1(j);
	}

	for (int j = 1; j <= N; j++)
	{
		fh(j, pt_test->gc)=0;
		fh(j, M + pt_test->gc)=0;
	}
	
	for (int i = pt_test->gc + 1; i <= M + pt_test->gc; i++)
		for (int j = 1; j <= N; j++)
    		K(j, i - pt_test->gc) = -(fh(j, i) - fh(j, i - 1))/h;

	delete ptA;
}