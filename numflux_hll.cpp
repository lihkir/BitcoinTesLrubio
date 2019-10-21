#include <math.h>
#include "fn_flux.h"
#include "numflux_hll.h"
#include "utilities.h"

std::vector<std::vector<double>> numflux_hll(std::vector<double>& ul, std::vector<double>& ur, std::vector<double>& Sl, std::vector<double>& Sr, std::vector<double>& S)
{
	int N = ul.size();

	S[0] = min(Sl[0], Sr[0]);
	S[1] = max(Sl[1], Sr[1]);

	std::vector<double> fl = fn_flux(ul);
	std::vector<double> fr = fn_flux(ur);

	std::vector<std::vector<double>> fh(N, std::vector<double>(1));

	if (S[0] >= 0) 
	{
		for (int i = 0; i < N; i++) 
			fh[i][0] = fl[i];
	}
	else if (S[1] <= 0)
	{
		for (int i = 0; i < N; i++) 
			fh[i][0] = fr[i];
	}
	else 
	{
		double ids = 1.0 / (S[1] - S[0]);
		double gp = S[1] * ids;
		double gm = -S[0] * ids;

		for (int i = 1; i < N; i++) 
			fh[i][0] = gp * (fl[i] - S[0] * ul[i]) + gm * (fr[i] - S[1] * ur[i]);
	}
	return fh;
}