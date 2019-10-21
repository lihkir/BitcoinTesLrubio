#include <math.h>
#include "bc.h"
#include "minmax_charspeed.h"
#include "muscl.h"
#include "numflux_hll.h"
#include "numflux_pvm.h"
#include "convec_pvm.h"
#include "utilities.h"

void convec_pvm(std::vector<std::vector<double>>& v, double h, std::vector<std::vector<double>>& K)
{
	struct test_cases* pt_test = get_tests();

	int M = v[0].size() - 2 *pt_test->gc; int N = v.size();

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
	
	std::vector<std::vector<double>> ul(N, std::vector<double>(M + 2 * pt_test->gc));
	std::vector<std::vector<double>> ur(N, std::vector<double>(M + 2 * pt_test->gc));
	std::vector<std::vector<double>> fh(N, std::vector<double>(M + 2 * pt_test->gc));

	std::vector<double> uli(N);
	std::vector<double> uri(N);
	std::vector<double> Sl(N);
	std::vector<double> Sr(N);

	auto N5 = max(N, 5);
	std::vector<std::vector<double>> A(N, std::vector<double>(N5));
	std::vector<std::vector<double>> Al(N, std::vector<double>(N));

	std::vector<double> fl(N);
	std::vector<double> fr(N);
	std::vector<double> ua(N);

	std::vector<std::vector<double>> dif(N, std::vector<double>(1));
	std::vector<std::vector<double>> Qurmul(N, std::vector<double>(1));
	std::vector<double> Qc(n_coef);
	std::vector<double> Si(2);

	bc(v);
	std::vector<std::vector<double>> S = minmax_charspeed(v);

	double cs = 0;
	for (int i = pt_test->gc; i < M + pt_test->gc; i++)
		cs = max(cs, max(fabs(S[0][i]), fabs(S[1][i])));

	muscl(v, ul, ur);

	for (int i = pt_test->gc; i < M + pt_test->gc - 1; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			uli[j] = ul[j][i];
			uri[j] = ur[j][i];
		}
		for (int j = 0; j < 2; j++) 
		{
			Sl[j] = S[j][i];
			Sr[j] = S[j][i + 1];
		}

		std::vector<std::vector<double>> fhn;
		if (global::idx_q > 0)
			fhn = numflux_pvm(uli, uri, Sl, Sr, global::idx_q, A, Al, Qurmul, Qc, Si);
		else
			fhn = numflux_hll(uli, uri, Sl, Sr, Si);
		for (int j = 0; j < N; j++)
			fh[j][i] = fhn[j][0];
	}

	for (int j = 0; j < N; j++)
	{
		fh[j][pt_test->gc - 1] = 0;
		fh[j][M + pt_test->gc - 1] = 0;
	}

	for (int i = pt_test->gc; i < M + pt_test->gc; i++)
		for (int j = 0; j < N; j++)
			K[j][i - pt_test->gc] = -(fh[j][i] - fh[j][i - 1]) / h;
}