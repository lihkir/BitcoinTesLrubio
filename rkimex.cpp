#include "apply_diffus.h"
#include "bc.h"
#include "charspeed.h"
#include "convec.h"
#include "diffus.h"
#include "lsolve.h"
#include "nwtsolve.h"
#include "rkimex.h"
#include "update_utilities.h"

void rkimex(Matrix<double> &A, Vector<double> &b, Matrix<double> &Ah, Vector<double> &bh, Matrix<double> &u0, double T, double cfl, double L, int convec_type0, int imex_type)
{
  /** convec_type=0 => PVM | convec_type=1 => GLF**/ 

  struct test_cases* pt_test = get_tests();
  
  int convec_type = convec_type0;

  int m = size(u0, 1);
  int n = size(u0, 2);

  Matrix<double> u(m, n + 2*pt_test->gc);
  update_inside(u, u0, pt_test->gc);
  
  bc(u);

  double h = L/n;
  double cs = charspeed(u);
  double dt = cfl*h/cs;

  double t = 0;
  int iter = 0;
  double next_t;
  
  while (t  < T)
  {
    next_t = t + dt;
    if (next_t > T)
    {
      dt = T - t;
      next_t = T;
    }

    if (imex_type == 0) {
      do_rkimex(A, b, Ah, bh, u, h, dt, 100, 1e-10);
    } else {
      do_lirkimex(A, b, Ah, u, h, dt);
    } 
    
    t = next_t;
    iter = iter + 1;
        
    if (imex_type == 0)
      printf("iter = %4d | t = %.16f | dt = %.16f | cs = %.16f\n", iter, t, dt, cs);
    else
      printf("iter = %4d | t = %.16f | dt = %.16f | cs = %.16f\n", iter, t, dt, cs);        
       
    cs = charspeed(u);
    dt = cfl*h/cs;
  }
  update_inside(u0, u, pt_test->gc);
}

void do_rkimex(Matrix<double> &A, Vector<double> &b, Matrix<double> &Ah, Vector<double> &bh, Matrix<double> &u, double h, double dt, int maxits, double tol)
{
  struct test_cases* pt_test = get_tests();

  int s = size(A, 2);
  int n = size(u, 2) - 2*pt_test->gc;
  int m = size(u, 1);

  Matrix<double> ul(m, n + 2*pt_test->gc);
  Matrix<double> K0(m, n);

  Block<double> K(s, m, n);
  Block<double> Kh(s + 1, m, n);

  convec(u, h, K0);

  for (int j = 1; j <= n; j++)
    for (int i = 1; i <= m; i++)
      Kh(1, i, j) = K0(i, j);

  Matrix<double> B(m, n);
  
  for (int l = 1; l <= s; l++)
  {
    for (int q = 1; q <= n; q++)
    {
      for (int p = 1; p <= m; p++)
      {
        B(p, q) = u(p, q + pt_test->gc);
        for (int j = 1; j <= l - 1; j++)
          B(p, q) = B(p, q) + dt*A(l, j)*K(j, p, q);
      }
    }
    
    for (int j = 1; j <= l; j++)
      for (int q = 1; q <= n; q++)
        for (int p = 1; p <= m; p++)
          B(p, q) = B(p, q) + dt*Ah(l + 1, j)*Kh(j, p, q);

    update_inside(ul, u, 0);
    nwtsolve(dt*A(l,l), B, ul, h, dt, maxits, tol);
    diffus(ul, h, K0);
  
    for (int j = 1; j <= n; j++)
      for (int i = 1; i <= m; i++)
        K(l, i, j) = K0(i, j);

    if (l < s || bh(s+1) != 0)
    {
      convec(ul, h, K0);
      for (int j = 1; j <= n; j++)
        for (int i = 1; i <= m; i++)
          Kh(l + 1, i, j) = K0(i, j);
    }
  }

  for (int q = 1; q <= n; q++)
    for (int p = 1; p <= m; p++)
      for (int j = 1; j <= s; j++)
        u(p, q + pt_test->gc) = u(p, q + pt_test->gc) + dt*(b(j)*K(j, p, q) + bh(j)*Kh(j, p, q));

  if (bh(s + 1) != 0)
  {
    for (int q = 1; q <= n; q++)
      for (int p = 1; p <= m; p++)
        u(p, q + pt_test->gc) = u(p, q + pt_test->gc) + dt*bh(s + 1)*Kh(s + 1, p, q);
  }

  bc(u);
}

/*****************************************************************************/

void do_lirkimex(Matrix<double> &A, Vector<double> &b, Matrix<double> &At, Matrix<double> &u, double h, double dt)
{
  struct test_cases* pt_test = get_tests();

  int s = size(A, 2);
  int n = size(u, 2) - 2*pt_test->gc;
  int m = size(u, 1);

  Matrix<double> K0(m, n);
  Block<double> K(s, m, n);

  Matrix<double> ut(m, n + 2*pt_test->gc);
  Matrix<double> uh(m, n + 2*pt_test->gc);
  Matrix<double> R(m, n);

  for (int l = 1; l <= s; l++)
  {
    for (int q = 1; q <= n; q++)
    {
      for (int p = 1; p <= m; p++)
      {
        uh(p, q + pt_test->gc) = u(p, q + pt_test->gc);
        ut(p, q + pt_test->gc) = u(p, q + pt_test->gc);
        for (int j = 1; j <= l - 1; j++)
        {
          uh(p, q + pt_test->gc) = uh(p, q + pt_test->gc) + dt*A (l, j)*K(j, p, q);
	        ut(p, q + pt_test->gc) = ut(p, q + pt_test->gc) + dt*At(l, j)*K(j, p, q);
        }
      }
    }
    bc(uh);
    bc(ut);
    /** K0=C(ut) **/ 
    convec(ut, h, K0);
    /** R=B(ut)*uh **/
    apply_diffus(ut, uh, h, R);
    
    for (int q = 1; q <= n; q++)
      for (int p = 1; p <= m; p++)
        K0(p, q) = K0(p, q) + R(p, q);

    /** K0=C(ut)+B(ut)*uh a la entrada **/
    lsolve(dt*A(l,l), ut, K0, h);
    /** K0=inv(I-dt*A(ll)*B(ut))*K0 a la salida **/

    for (int j = 1; j <= n; j++)
      for (int i = 1; i <= m; i++)
        K(l, i, j)  = K0(i, j);
  }

  for (int q = 1; q <= n; q++)
    for (int p = 1; p <= m; p++)
      for (int j = 1; j <= s; j++)
        u(p, q + pt_test->gc) = u(p, q + pt_test->gc) + dt*b(j)*K(j, p, q);

  bc(u);
}