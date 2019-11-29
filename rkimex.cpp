#include <stdio.h>
#include <math.h>
#include "apply_diffus.h"
#include "bc.h"
#include "charspeed.h"
#include "convec.h"
#include "diffus.h"
#include "lsolve.h"
#include "nwtsolve.h"
#include "rkimex.h"
#include "utilities.h"
#include "save_data.h"

void rkimex(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& Ah, std::vector<double>& bh, std::vector<std::vector<double>>& u0, std::vector<double> Ta, double cfl, double L, int convec_type, int imex_type, int test)
{
  struct test_cases* pt_test = get_tests();
  double relres; 
  
  int m = u0.size();
  int n = u0[0].size();

  std::vector<std::vector<double>> u(m, std::vector<double>(n + 2*pt_test->gc));
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      u[i][pt_test->gc + j] = u0[i][j];
  
  bc(u);

  double h = L/n;
  double cs = charspeed(u);

  double dt = test == 2 ? 500*h : cfl*h/cs;

  double t = 0, T = 0;
  int iter = 0;
  double next_t;
  
  clock_t start, end;
  start = clock();
  for (unsigned int i = 0; i < Ta.size(); i++) 
  {
    T = Ta[i];
    while (t <  T) 
    {
      next_t = t + dt;
      if  (next_t > T) 
      {
        dt = T - t;
        next_t = T;
      }

      if  (imex_type == 0)
        do_rkimex(A, b, Ah, bh, u, h, dt, 100, 1e-10);
      else
        relres = do_lirkimex(A, b, Ah, u, h, dt);
        
      t = next_t;
      iter = iter + 1;
      printf("iter:%4d => t: %0.14f\t dt: %0.14f\t cs: %0.14f\t relres = %.14f\n", iter, t, dt, cs, relres);
    
      cs=charspeed(u);
      dt = test == 2 ? 500*h : cfl*h/cs;
    }
    end = clock();
  	double CPUTIME = double(end - start) / CLOCKS_PER_SEC;
    printf("\n##########################################\n");
  	printf("\nCPUTIME rkimex T=%.15f=> %.15f\n", T, CPUTIME);
    printf("\n##########################################\n");
    
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        u0[i][j] = u[i][pt_test->gc + j];
    
    save_data(u0, pt_test->gc, T);
  }
}

void do_rkimex(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& Ah, std::vector<double>& bh, std::vector<std::vector<double>>& u, double h, double dt, int maxits, double tol)
{
  struct test_cases* pt_test = get_tests();

  int s = A[0].size();
  int n = u[0].size() - 2 * pt_test->gc;
  int m = u.size();

  std::vector<std::vector<double>> ul(m, std::vector<double>(n + 2 * pt_test->gc));
  std::vector<std::vector<double>> K0(m, std::vector<double>(n));
  std::vector<std::vector<double>> Mele(m , std::vector<double>(n));

  std::map<int, std::vector<std::vector<double>>> K;
	for (int i = 0; i < s; i++) K[i] = Mele;

  std::map<int, std::vector<std::vector<double>>> Kh; 
  for (int i = 0; i < s+1; i++) Kh[i] = Mele;

  convec(u, h, K0);
  for (int j = 0; j < n; j++) 
    for (int i = 0; i < m; i++ )
      Kh[0][i][j] = K0[i][j];

  std::vector<std::vector<double>> B(m, std::vector<double>(n));

  for (int l = 0; l < s; l++)
  {
    for (int q = 0; q < n; q++)
    {
      for (int p = 0; p < m; p++)
      {
        B[p][q] = u[p][q+pt_test->gc];
        for (int j = 0; j < l-1; j++) B[p][q] = B[p][q] + dt*A[l][j]*K[j][p][q];
       }
    }
    
    for (int j = 0; j < l; j++)
      for (int q = 0; q < n; q++)
        for (int p = 0; p < m; p++)
          B[p][q] = B[p][q] + dt*Ah[l+1][j]*Kh[j][p][q];

    std::vector<std::vector<double>> ul = SubMatrix(u, 0);  
    nwtsolve(dt*A[l][l], B, ul, h, dt, maxits, tol); 
    diffus(ul, h, K0);

    for (int j = 0; j < n; j++) 
      for (int i = 0; i < m; i++)
        K[l][i][j] = K0[i][j];

    if  (l < s || bh[s + 1] !=0)
    {
      convec(ul, h, K0);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          Kh[l + 1][i][j] = K0[i][j];
    }
  }

  for (int q = 0; q < n; q++) 
    for (int p = 0; p < m; p++) 
      for (int j = 0; j < s; j++) 
        u[p][q + pt_test->gc] = u[p][q + pt_test->gc] + dt*(b[j]*K[j][p][q] + bh[j]*Kh[j][p][q]);

  if (bh[s+1] !=0)
  {
    for (int q = 0; q < n; q++) 
      for (int p = 0; p < m; p++) 
        u[p][q + pt_test->gc] = u[p][q + pt_test->gc] + dt*bh[s+1]*Kh[s+1][p][q];
  }
  bc(u);
}

double do_lirkimex(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& At, std::vector<std::vector<double>>& u, double h, double dt)
{
  struct test_cases* pt_test = get_tests();
  
  double relres = 0; 
  int s = A[0].size();
  int n = u[0].size() - 2*pt_test->gc;
  int m = u.size();

  std::map<int, std::vector<std::vector<double>>> K;
  std::vector<std::vector<double>> Mele(m , std::vector<double>(n));
	for (int i = 0; i < s; i++) K[i] = Mele;

  std::vector<std::vector<double>> ut(m, std::vector<double>(n + 2*pt_test->gc));
  std::vector<std::vector<double>> uh(m, std::vector<double>(n + 2*pt_test->gc));  
  std::vector<std::vector<double>> K0(m, std::vector<double>(n));
  std::vector<std::vector<double>>  R(m, std::vector<double>(n));

  for (int l = 0; l < s; l++)
  {
    for (int q = 0; q < n; q++) 
    {
      for (int p = 0; p < m; p++) 
      {
        uh[p][q + pt_test->gc] = u[p][q + pt_test->gc];
        ut[p][q + pt_test->gc] = u[p][q + pt_test->gc];
        for (int j = 0; j < l-1; j++) 
        {
          uh[p][q + pt_test->gc] = uh[p][q + pt_test->gc] + dt*A[l][j]*K[j][p][q];
          ut[p][q + pt_test->gc] = ut[p][q + pt_test->gc] + dt*At[l][j]*K[j][p][q];
        }
      }
    }
    bc(uh);
    bc(ut);
    ////////////////////// K0=C(ut)
    convec(ut, h, K0);
    ////////////////////// R=B(ut)*uh
    apply_diffus(ut, uh, h, R);
    for (int q = 0; q < n; q++)
      for (int p = 0; p < m; p++)
        K0[p][q] = K0[p][q] + R[p][q];
    ////////////////// K0=C(ut)+B(ut)*uh a la entrada
    relres = lsolve(dt*A[l][l], ut, K0, h);
    ////////////////// K0=inv(I-dt*A(ll)*B(ut))*K0 a la salida
    for (int j = 0; j < n; j++) 
      for (int i = 0; i < m; i++)
        K[l][i][j] = K0[i][j];
  }

  for (int q = 0; q < n; q++)
    for (int p = 0; p < m; p++)
      for (int j = 0; j < s; j++)
        u[p][q + pt_test->gc] = u[p][q + pt_test->gc] + dt*b[j]*K[j][p][q];

  bc(u);
  return relres;
}