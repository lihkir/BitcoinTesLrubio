#include <math.h>
#include "fn_flux.h"
#include "numflux_hll.h"
#include "update_utilities.h"

void numflux_hll(Matrix<double> &ul, Matrix<double> &ur, Matrix<double> &Sl, Matrix<double> &Sr, Matrix<double> &fh, Matrix<double> &fl, Matrix<double> &fr, Matrix<double> &S)
{
	int N = size(ul, 1);
	double ids, gp, gm;

	S(1) = min(Sl(1), Sr(1));
	S(2) = max(Sl(2), Sr(2));

	fn_flux(ul,fl);
	fn_flux(ur,fr);

	if (S(1) >= 0)
		update_inside(fh, fl, 0);
	else if (S(2) <= 0)
		update_inside(fh, fr, 0);
 	else
	{
		ids = 1/(S(2) - S(1));
  		gp = S(2)*ids;
  		gm = -S(1)*ids;

		for (int i = 1; i <= N; i++)
			fh(i) = gp*(fl(i) - S(1)*ul(i)) + gm*(fr(i) - S(2)*ur(i));
	}
}