#include <stdio.h>
#include <math.h>
#include "residual.h"

void residual(Block<double> &D, Block<double> &L, Block<double> &U, Matrix<double> &v, Matrix<double> &r)
{
  int m = size(r, 1);
  int n = size(D, 1);

  for (int i = 1; i <= m; i++)
  {
    for (int k = 1; k <= m; k++)
      r(i, 1) = r(i, 1) - D(1, i, k)*v(k, 1) - U(1, i, k)*v(k, 2);

    for (int j = 2; j <= n - 1; j++)
      for (int k = 1; k <= m; k++)
        r(i, j) = r(i, j) - L(j-1, i, k)*v(k,j - 1) - D(j, i, k)*v(k, j) - U(j, i, k)*v(k,j + 1);
    
    for (int k = 1; k <= m; k++)
      r(i, n) = r(i, n) - L(n - 1, i, k)*v(k, n-1) - D(n, i, k)*v(k, n);
  }
}