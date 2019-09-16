#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <pthread.h>
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "average.h"
#include "utilities.h"
#include "bc.h"
#include "fn_flux.h"
#include "matmult.h"
#include "hornerm.h"
#include "test_cases.h"
#include "lminmax_charspeed.h"

using namespace std;

extern pthread_mutex_t lock;
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; 

int main(int argc, char* argv[])
{
  // M_rows := number of species and N_cols := number of nodes
  int N_cols  = atoi(argv[1]);

  extern int idx_test;
  pthread_mutex_lock(&lock);
  idx_test = atoi(argv[2]);
  pthread_mutex_unlock(&lock);

  srand (time(NULL));

  printf("\n#######################################################################\n");
  printf("\ntesting: beta_cases.cpp function\n\n");

  struct test_cases *pt_test = get_tests();
  printf("pt_test->M_rows = %d\n\n", pt_test->M_rows);
  PrintingContainer(pt_test->beta);

  printf("\n#######################################################################\n");
  printf("\ntesting: average.cpp function\n\n");

  std::vector<double> ul = RandomVector<double>(pt_test->M_rows, 9);
  std::vector<double> ur = RandomVector<double>(pt_test->M_rows, 9);
  std::vector<double> az = RandomVector<double>(pt_test->M_rows, 9);

  PrintingContainer(ul);
  PrintingContainer(ur);
  std::vector<double> ua = average(ul, ur, 0.5);
  PrintingContainer(ua);

  printf("\n#######################################################################\n");
  printf("\nTesting: bc.cpp function\n\n");

  std::vector<std::vector<double>> m   = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);
  std::vector<std::vector<double>> z   = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);
  std::vector<std::vector<double>> p   = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);
  std::vector<std::vector<double>> A   = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);
  std::vector<std::vector<double>> x   = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);
  std::vector<std::vector<double>> pAx = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);
  std::vector<std::vector<double>> pr  = RandomMatrix<double>(pt_test->M_rows, N_cols, 9);

  printf("\nm, z, p matrices before ghost:\n\n");
  PrintingContainer(m);
  PrintingContainer(z);
  PrintingContainer(p);
  bc(m);
  bc(z);
  bc(p);
  printf("\nm, z, p matrices after ghost:\n\n");
  PrintingContainer(m);
  PrintingContainer(z);
  PrintingContainer(p);

  printf("testing fn_flux.cpp function\n\n");

  std::vector<double> fn = fn_flux(ul, pt_test->beta);
  PrintingContainer(fn);
  
  printf("#######################################################################\n\n");
  printf("\ntesting hornerm.cpp\n");

  bc(A);
  bc(x);
  bc(pr);
  bc(pAx);
  printf("\naz, A, x, pr, pAx: before\n");
  PrintingContainer(az);
  PrintingContainer(A);
  PrintingContainer(x);
  PrintingContainer(pr);
  PrintingContainer(pAx);
  hornerm(A, x, az, pr, pAx);
  printf("\naz, A, x, pr, pAx: later\n");
  PrintingContainer(az);
  PrintingContainer(A);
  PrintingContainer(x);
  PrintingContainer(pr);
  PrintingContainer(pAx);
  
  printf("#######################################################################\n\n");  
  
  double sed_hsfd = sed_hsf(0.5);
  printf("testing: sed_hsfd = %f\n\n", sed_hsfd);
  
  double sed_hsfd_der = sed_hsf_der(1);
  printf("testing: sed_hsfd = %f\n\n", sed_hsfd_der);
  
  printf("#######################################################################\n\n");
  printf("testing matmult.cpp function\n\n");

  printf("\nm, z, p before:\n");
  PrintingContainer(m);
  PrintingContainer(z);
  PrintingContainer(p);
  matmult(m, z, p);
  printf("\nm, z, p later:\n");
  PrintingContainer(m);
  PrintingContainer(z);
  PrintingContainer(p);

  printf("#######################################################################\n\n");
  printf("testing lminmax_charspeed.cpp function\n\n");

  std::vector<double> S = lminmax_charspeed(ul, ur);

  PrintingContainer(S);

  return 0;
}