#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iterator>  
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "utilities.h"

std::vector<double> lminmax_charspeed(std::vector<double>& phil, std::vector<double>& phir)
{
	struct test_cases* pt_test = get_tests();
	auto N = pt_test->delta.size();

	std::vector<double> S;
	double s_phil = 0, s_phir = 0, etal = 0, etar = 0;

	for (auto j = 0; j < N; j++)
	{
		s_phil += phil[j];
		etal += pt_test->delta[j] * phil[j];
		s_phir += phir[j];
		etar += pt_test->delta[j] * phir[j];
	}

	double s_delta = s_phil - s_phir;
	double eta_delta = etal - etar;

	double w1 = sed_hsf(s_phil);
	double w2 = sed_hsf(s_phir);

	double maxv1 = pt_test->delta[0] * w1;
	double maxv2 = pt_test->delta[0] * w2;
	double max_0 = max(maxv1, maxv2);

	std::vector<double> buf_opt;
	buf_opt.push_back(w1 * pt_test->delta[N - 1] + etal * sed_hsf_der(s_phil));
	buf_opt.push_back(w2 * pt_test->delta[N - 1] + etar * sed_hsf_der(s_phir));

	if (eta_delta != 0 && s_delta != 0)
	{
		double Q = pt_test->nexp * eta_delta * s_delta + ((s_delta) * (s_delta)) * pt_test->delta[N - 1];
		double alpha = ((1 - s_phir) * (pt_test->delta[N - 1] * s_delta + eta_delta) - (pt_test->nexp - 1) * etar * s_delta) / Q;
		if (alpha > 0 && alpha < 1)
		{
			double vphi = s_phir + alpha * s_delta;
			double etav = etar + alpha * eta_delta;
			double w = sed_hsf(vphi);
			buf_opt.push_back(w * pt_test->delta[N - 1] + etav * sed_hsf_der(vphi));
		}
	}

	auto min_n = *min_element(begin(buf_opt), end(buf_opt));
	S.push_back(min_n);
	S.push_back(max_0);

	return S;
}