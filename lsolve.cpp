#include <stdio.h>
#include "diffusion_matrix.h"
#include "residual.h"
#include "trdsolve.h"
#include "lsolve.h"
#include "test_cases.h"
#include "update_utilities.h"

double lsolve(double a, Matrix<double> &u, Matrix<double> &b, double h)
{
  struct test_cases* pt_test = get_tests();

  int m = size(u, 1);
  int n = size(u, 2) - 2*pt_test->gc;

  Matrix<double> v(m, n);
  Matrix<double> bj(m, 1);
  Matrix<double> rj(m, 1);
  update_inside(v, u, pt_test->gc);
  
  Block<double> D(n, m, m);
  Block<double> L(n, m, m);
  Block<double> U(n, m, m);

  diffusion_matrix(v, h, a, D, L, U);
  
  double nrm2_b = 0;
  for (int j = 1; j <= n; j++)
  {
    get_col(bj, b, j);
    nrm2_b += dot_product(bj, bj);
  }
  
  Matrix<double> r(m, n);

  update_inside(v, b, 0);
  trdsolve(D, L, U, b);

  update_inside(r, v, 0);
  residual(D, L, U, b, r);

  double nrm2_r=0;
  for (int j = 1; j <= n; j++)
  {
    get_col(rj, r, j);
    nrm2_r +=  dot_product(rj, rj);
  }

  return sqrt(nrm2_r/nrm2_b);
}