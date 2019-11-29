#include <math.h>
#include "average.h"
#include "fn_flux.h"
#include "hornerm.h"
#include "jacobiana.h"
#include "jacobiana_dec.h"
#include "numflux_pvm.h"
#include "Qcoeff.h"
#include "utilities.h"

std::vector<std::vector<double>> numflux_pvm(std::vector<double>& ul, std::vector<double>& ur, std::vector<double>& Sl, std::vector<double>& Sr, int idx_q, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Al, std::vector<std::vector<double>>& Qurmul, std::vector<double>& Qc, std::vector<double>& S)
{
	int N = ul.size();
	
	std::vector<std::vector<double>> dif(N, std::vector<double>(1));
	for (int i = 0; i < N; i++) dif[i][0] = ur[i] - ul[i];

	if (global::int_form == 1) 
	{
		/** Midpoint Rule **/
		std::vector<double> ua = average(ul, ur, 0.5);
		std::vector<std::vector<double>> A = jacobiana(ua);
	}
	else if (global::int_form == 5) 
	{
		/** Midpoint Rule deconstructed Jacobian **/
		std::vector<double> ua = average(ul, ur, 0.5);
		std::vector<std::vector<double>> A = jacobiana_dec(ua);
	}
	else if (global::int_form == 2) 
	{
		/** Gaussian Rule 2 nodes **/
		/** 0.211324865405187=0.5*(1-1.0/sqrt(3)) **/
		std::vector<double> ua1 = average(ul, ur, 0.211324865405187);
		std::vector<std::vector<double>> Al1 = jacobiana(ua1);

		for (int i = 0; i < N; i++) 
			for (int j = 0; j < N; j++) 
				A[i][j] = 0.5 * Al1[i][j];

		/** 0.788675134594813=0.5*(1+1.0/sqrt(3)) **/
		std::vector<double> ua2 = average(ul, ur, 0.788675134594813);
		std::vector<std::vector<double>> Al2 = jacobiana(ua2);

		for (int i = 0; i < N; i++) 
			for (int j = 0; j < N; j++) 
				A[i][j] = A[i][j] + 0.5 * Al2[i][j];
	}
	else if (global::int_form == 3) 
	{
		/** Gaussian Rule 3 nodes **/
		std::vector<double> ua1 = average(ul, ur, 0.112701665379258);
		std::vector<std::vector<double>> Al1 = jacobiana(ua1);

		for (int i = 0; i < N; i++) 
			for (int j = 0; j < N; j++) 
				A[i][j] = 5.0 / 18 * Al1[i][j];

		std::vector<double> ua2 = average(ul, ur, 0.112701665379258);
		std::vector<std::vector<double>> Al2 = jacobiana(ua2);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] + 5.0 / 18 * Al2[i][j];
			
		std::vector<double> ua3 = average(ul, ur, 0.5);
		std::vector<std::vector<double>> Al3 = jacobiana(ua3);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) 
				A[i][j] = A[i][j] + 8.0 / 18 * Al3[i][j];
	}
	else if (global::int_form == 4)
	{
		std::vector<std::vector<double>> A(N, std::vector<double>(N));
		double ds = 1e-4;
		for (double s = ds / 2; s <= 1; s += ds) {
			
			std::vector<double> ua = average(ul, ur, s);
			std::vector<std::vector<double>> Al = jacobiana(ua);

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					A[i][j] = A[i][j] + Al[i][j];
		}
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				A[i][j] = ds * A[i][j];
	}

	S[0] = min(Sl[0], Sr[0]);
	S[1] = max(Sl[1], Sr[1]);

	std::vector<double> fl = fn_flux(ul);
	std::vector<double> fr = fn_flux(ur);
	std::vector<std::vector<double>> fh(N, std::vector<double>(1));

	for (unsigned int j = 0; j < fh[0].size(); j++)
	{
		if (S[0] >= 0) 
		{
			for (int i = 0; i < N; i++)
				fh[i][j] = fl[i];
		}
		else if (S[1] <= 0) 
		{
			for (int i = 0; i < N; i++)
				fh[i][j] = fr[i];
		}
		else
		{
			std::vector<double> Qc = Qcoeff(S, idx_q);
			hornerm(A, dif, Qc, fh, Qurmul);
			for (int i = 0; i < N; i++)
				fh[i][j] = 0.5 * (fl[i] + fr[i] - Qurmul[i][j]);
		}
	}
	return fh;
}