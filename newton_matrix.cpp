#include <stdio.h>
#include <math.h>
#include "der_diffusion_tensor.h"
#include "diffusion_tensor.h"
#include "newton_matrix.h"

void newton_matrix(Matrix<double> &u, double h, double a, Block<double> &D, Block<double> &L, Block<double> &U)
{
    int n = size(u, 2);
    int m = size(u, 1);

    double f = -a*0.5/(h*h);
    double z;    

    Block<double> B(n, m, m);
    Block<double> Bp(n, m, m);
    diffusion_tensor(u, B);

    for (int j = 1; j <= n; j++)
    {
        for (int i = 1; i <= m; i++)
        {
            for (int p = 1; p <= m; p++)
            {
                D(j, p, i)=0;
                L(j, p, i)=0;
                U(j, p, i)=0;
            }
        }
    }

    for (int r = 1; r <= m; r++)
    {
        der_diffusion_tensor(u, Bp, r);
        for (int p = 1; p <= m; p++)
        {
            for (int i = 2; i <= n; i++)
            {
                z = B(i, p, r) + B(i -1, p, r);
                D(i, p, r) = D(i, p, r) - f*z;
                L(i - 1, p, r) = L(i - 1, p, r) + f*z;
                
                for (int q = 1; q <= m; q++)
                    D(i, p, r) = D(i, p, r) - f*Bp(i, p, q)*(u(q, i)-u(q, i-1));

                for (int q = 1; q <= m; q++)
                    L(i - 1, p, r) = L(i - 1, p, r) - f*Bp(i - 1, p, q)*(u(q, i)-u(q, i-1));
            }
            for (int i = 1; i <= n - 1; i++)
            {
                z = B(i + 1, p, r) + B(i, p, r);
                U(i, p, r) = U(i, p, r) + f*z;
                D(i, p, r) = D(i, p, r) - f*z;
                
                for (int q = 1; q <= m; q++)
                    D(i, p, r) = D(i, p, r) + f*Bp(i, p, q)*(u(q, i+1)-u(q, i));
                
                for (int q = 1; q <= m; q++)
                    U(i, p, r) = U(i, p, r) + f*Bp(i + 1, p, q)*(u(q, i+1)-u(q, i));
            }
        }
    }

    for (int i = 1; i <= n; i++)
        for (int r = 1; r <= m; r++)
            D(i, r, r) = D(i, r, r) + 1;
}