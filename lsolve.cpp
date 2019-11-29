#include <stdio.h>
#include <math.h>
#include "diffusion_matrix.h"
#include "residual.h"
#include "trdsolve.h"
#include "lsolve.h"
#include "test_cases.h"
#include "utilities.h"

double lsolve(double a, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &b, double h)
{
  struct test_cases* pt_test = get_tests();

  int m = u.size();
  int n = u[0].size() - 2 * pt_test->gc;

  std::vector<std::vector<double>> v = SubMatrix(u, pt_test->gc);
  
  std::map<int, std::vector<std::vector<double>>> D;
  std::map<int, std::vector<std::vector<double>>> L;
  std::map<int, std::vector<std::vector<double>>> U;
  
  diffusion_matrix(v, h, a, D, L, U);

  double nrm2_b=0;
  for (int j = 0; j < n; j++) nrm2_b += SquareSum(SubVector(b, j));

  std::vector<std::vector<double>> w = SubMatrix(b, 0);
  trdsolve(D, L, U, b);
  std::vector<std::vector<double>> r = SubMatrix(w, 0);
  residual(D, L, U, b, r);
  
  double nrm2_r = 0;
  for (int j = 0; j < n; j++) nrm2_r += SquareSum(SubVector(r, j));

  return sqrt(nrm2_r/nrm2_b);
}