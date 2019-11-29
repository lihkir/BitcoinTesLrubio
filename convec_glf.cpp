#include <stdio.h>
#include <math.h>
#include "charspeed.h"
#include "fn_flux.h"
#include "weno5.h"
#include "convec_glf.h"
#include "test_cases.h"
#include "utilities.h"

void convec_glf(std::vector<std::vector<double>>& v, double h, std::vector<std::vector<double>>& K)
{
	struct test_cases *pt_test = get_tests();

	int m = v.size();
	int n = v[0].size() - 2 * pt_test->gc;

	std::vector<std::vector<double>> fh(m , std::vector<double>(n+1));
	numflux_glf(v, fh);
	
	for (int i = 0; i < m; i++) 
		for (int j = 0; j < n; j++) 
			K[i][j] = -(fh[i][j+1] - fh[i][j]) / h; /** fh(*, j)=> fh en j+2.5=> K(*, j)=Dfh en j+3 **/ 
}

void numflux_glf(std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& fh)
{
	struct test_cases *pt_test = get_tests();

	int m = v.size();
	int n = v[0].size() - 2 * pt_test->gc;

	double cs = charspeed(v);

	std::vector<double> fi(m);
	std::vector<double> vi(m);
	std::vector<std::vector<double>> f(m, std::vector<double>(n + 2*pt_test->gc));

	for (int i = 0; i < n + 2 * pt_test->gc; i++)
	{
		vi = SubVector(v, i); fi = fn_flux(vi);
		for (int j = 0; j < m; j++) 
			f[j][i] = fi[j];
	}

	std::vector<std::vector<double>> fp(m, std::vector<double>(n + 2*pt_test->gc));
	std::vector<std::vector<double>> fm(m, std::vector<double>(n + 2*pt_test->gc));

	for (int i = 0; i < m; i++) 
	{
		for (int j = 0; j < n + 2 * pt_test->gc; j++) 
		{
			fp[i][j] = 0.5 * (f[i][j] + cs * v[i][j]);
			fm[i][j] = 0.5 * (f[i][j] - cs * v[i][j]);
		}
	}

	std::vector<double> f1(5);
	for (int k = 0; k < m; k++)
	{
		for (int i = 1; i < n; i++) 
		{
			for (int j = 0; j < 5; j++) f1[j] = fp[k][i+j-1]; /** Stencil: i+4 flujo+ en i+2.5 **/
			double fp0 = weno5(f1, 1e-100);

			for (int j = 0; j < 5; j++) f1[j] = fm[k][i + 6 - j]; /** Flujo en i+2.5 **/ 
			double fm0 = weno5(f1, 1e-100);

			fh[k][i] = fp0 + fm0; /** fh(*, i)=> Flujo en i+2.5 **/
		}
	}
	
	for (int k = 0; k < m; k++) {
		fh[k][0] = 0; fh[k][n] = 0;
	}
}