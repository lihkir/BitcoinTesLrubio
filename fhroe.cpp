#include "fhroe.h"
#include "fn_flux.h"

void fhroe(Matrix<double> &ul, Matrix<double> &ur, Matrix<double> &fl, Matrix<double> &fr, Matrix<double> &fh)
{
    fn_flux(ul,fl);
    fn_flux(ur,fr);

    int N = size(fh, 1);
    for (int i = 1; i <= N; i++)
        fh(i) = 0.5*(fl(i) + fr(i) - abs(fr(i) - fl(i))*sgn(ur(i) - ul(i)));
}