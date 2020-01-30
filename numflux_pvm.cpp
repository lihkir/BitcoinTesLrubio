#include <math.h>
#include "average.h"
#include "fn_flux.h"
#include "hornerm.h"
#include "jacobiana.h"
#include "jacobiana_dec.h"
#include "numflux_pvm.h"
#include "Qcoeff.h"
#include "update_utilities.h"

void numflux_pvm(Matrix<double> &ul, Matrix<double> &ur, Matrix<double> &Sl, Matrix<double> &Sr, int idx_q, Matrix<double> &fh, Matrix<double> &A, Matrix<double> &Al, Matrix<double> &fl, Matrix<double> &fr, Matrix<double> &ua, Matrix<double> &dif, Matrix<double> &Qurmul, Matrix<double> &Qc, Matrix<double> &S)
{
	int N = size(ul, 1);
	double ds;

	for (int i = 1; i <= N; i++)
		dif(i) = ur(i) - ul(i);

	if (global::int_form == 1)
	{/** Midpoint Rule **/
		average(ul, ur, 0.5, ua);
  		jacobiana(ua, A);
	}
	else if (global::int_form == 5)
	{/** Midpoint Rule deconstructed Jacobian **/
  		average(ul, ur, 0.5, ua);
  		jacobiana_dec(ua, A);
	}
	else if (global::int_form == 2)
	{/** Gaussian Rule 2 nodes: 0.211324865405187=0.5*(1-1/sqrt(3)) **/
  		average(ul, ur, 0.211324865405187, ua);
  		jacobiana(ua, Al);
		
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
				A(i,j)=0.5*Al(i,j);
	
		/** 0.788675134594813=0.5*(1+1/sqrt(3)) **/
  		average(ul, ur, 0.788675134594813, ua);
  		jacobiana(ua, Al);
  
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
				A(i,j)=A(i, j)+0.5*Al(i,j);
	}  
	else if (global::int_form == 3)
	{/** Gaussian Rule 3 nodes **/
  		average(ul, ur, 0.112701665379258, ua);
  		jacobiana(ua, Al);

		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
      			A(i, j) = 5/18*Al(i,j);

  		average(ul, ur,  0.887298334620742, ua);
  		jacobiana(ua, Al);
  
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
      			A(i, j) = A(i, j) + 5/18*Al(i,j);

  		average(ul, ur, 0.5, ua);
  		jacobiana(ua, Al);

		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
    			A(i,j) = A(i, j) + 8/18*Al(i,j);
	}
	else if (global::int_form == 4)
	{
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
      			A(i, j) = 0;
		ds = 1e-4;
		for (double s = ds/2; s <= 1; s += ds)
		{
			average(ul, ur, s, ua);
    		jacobiana(ua, Al);

			for (int i = 1; i <= N; i++)
				for (int j = 1; j <= N; j++)
					A(i, j) = A(i, j) + Al(i,j);
		}      
		
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
				A(i, j) = ds*A(i, j);
	}

	S(1) = min(Sl(1), Sr(1));
	S(2) = max(Sl(2), Sr(2));
	
	fn_flux(ul, fl);
	fn_flux(ur, fr);

	if (S(1) >= 0)
		update_inside(fh, fl, 0);
	else if (S(2) <= 0)
		update_inside(fh, fr, 0);
  	else
	  {
		  Qcoeff(S, idx_q, Qc);
		  hornerm(A, dif, Qc, fh, Qurmul);
		  for (int i = 1; i <= N; i++)
		  	fh(i) = 0.5*(fl(i) + fr(i) - Qurmul(i));
	  }
}