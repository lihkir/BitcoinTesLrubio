#include <stdio.h>
#include <math.h>
#include "charspeed.h"
#include "fn_flux.h"
#include "weno5.h"
#include "convec_glf.h"
#include "test_cases.h"
#include "update_utilities.h"

void convec_glf(Matrix<double>& v, double h, Matrix<double>& K)
{
	struct test_cases* pt_test = get_tests();
	
	int n = size(v, 2) - 2*pt_test->gc;
	int m = size(v, 1);

	Matrix<double> *ptfh = new Matrix<double>(m, n + 1); Matrix<double> &fh = *ptfh;
	numflux_glf(v, fh);

	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{/** fh(*, j) => fh en j+2.5 => K(*, j)= Dfh en j+3 **/
			K(i, j) = -(fh(i, j + 1) - fh(i, j))/h;
		}
	}	
	delete ptfh;
}

void numflux_glf(Matrix<double>& v, Matrix<double>& fh)
{
	struct test_cases* pt_test = get_tests();

	double cs = charspeed(v);
	int n = size(v, 2) - 2*pt_test->gc;
	int m = size(v, 1);
	double fp0, fm0;

	Matrix<double> *ptfi = new Matrix<double>(m, 1); Matrix<double> &fi = *ptfi;
	Matrix<double> *ptvi = new Matrix<double>(m, 1); Matrix<double> &vi = *ptvi;
	Matrix<double> *ptf  = new Matrix<double>(m, n + 2*pt_test->gc); Matrix<double> &f = *ptf;
	
	for (int i = 1; i <= n + 2*pt_test->gc; i++)
	{
		get_col(vi, v, i);
  		fn_flux(vi, fi);
		for (int j = 1; i <= m; i++)
			f(j, i) = fi(j);
	}

	Matrix<double> *ptfp = new Matrix<double>(m, n + 2*pt_test->gc); Matrix<double> &fp = *ptfp;
	Matrix<double> *ptfm = new Matrix<double>(m, n + 2*pt_test->gc); Matrix<double> &fm = *ptfm;

	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n + 2*pt_test->gc; j++)
		{
			fp(i, j) = 0.5*(f(i, j) + cs*v(i, j));
    		fm(i, j) = 0.5*(f(i, j) - cs*v(i, j));
		}
	}

	Matrix<double> *ptf1 = new Matrix<double>(5, 1); Matrix<double> &f1 = *ptf1;
	
	for (int k = 1; k <= m; k++)
	{
		for (int i = 2; i <= n; i++)
		{
			for (int j = 1; j <= 5; j++)
      			f1(j) = fp(k, i+j-1); 
	
    		fp0 = weno5(f1, 1e-100);
	
			for (int j = 1; j <= 5; j++)
				f1(j) = fm(k, i+6-j); 

    		fm0 = weno5(f1, 1e-100);
			fh(k, i) = fp0+fm0;
		}
	}

	for (int i = 1; i <= m; i++)
	{
		fh(i, n + 1) = 0;
		fh(i, 1) = 0;
	}

	delete ptf;
	delete ptf1;
	delete ptfi;
	delete ptfm;
	delete ptfp;
	delete ptvi;
}