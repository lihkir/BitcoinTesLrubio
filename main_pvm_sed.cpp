#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "globals.h"
#include "sed_hsf.h"
#include "average.h"
#include "utilities.h"
#include "beta_cases.h"
#include "bc.h"
#include "fn_flux.h"
#include "matmult.h"
#include "hornerm.h"

std::vector<double>* pt_beta;

using namespace std;

int main(int argc, char* argv[])
{
  /** M_rows := number of species and N_cols := number of nodes **/

  int N_cols  = atoi(argv[1]);
  int test_id = atoi(argv[2]);

  global_data *pt_data;
  srand (time(NULL));

  printf("\n#######################################################################\n");
  printf("\ntesting: beta_cases.cpp function\n\n");

  int M_rows = 0;
  std::vector<double> beta = beta_cases(test_id, M_rows, pt_beta);
  PrintingContainer(beta);

  printf("\n#######################################################################\n");
  printf("\ntesting: average.cpp function\n\n");

  std::vector<double> ul = RandomVector<double>(M_rows, 9);
  std::vector<double> ur = RandomVector<double>(M_rows, 9);
  std::vector<double> az = RandomVector<double>(M_rows, 9);

  PrintingContainer(ul);
  PrintingContainer(ur);
  std::vector<double> ua = average(ul, ur, 0.5);
  PrintingContainer(ua);

  printf("\n#######################################################################\n");
  printf("\nTesting: bc.cpp function\n\n");

  std::vector<std::vector<double>> m   = RandomMatrix<double>(M_rows, N_cols, 9);
  std::vector<std::vector<double>> z   = RandomMatrix<double>(M_rows, N_cols, 9);
  std::vector<std::vector<double>> p   = RandomMatrix<double>(M_rows, N_cols, 9);
  std::vector<std::vector<double>> A   = RandomMatrix<double>(M_rows, N_cols, 9);
  std::vector<std::vector<double>> x   = RandomMatrix<double>(M_rows, N_cols, 9);
  std::vector<std::vector<double>> pAx = RandomMatrix<double>(M_rows, N_cols, 9);
  std::vector<std::vector<double>> pr  = RandomMatrix<double>(M_rows, N_cols, 9);

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

  printf("#######################################################################\n\n");  

  double sed_hsfd = sed_hsf(0.5);
  printf("testing: sed_hsfd = %f\n\n", sed_hsfd);
  
  printf("#######################################################################\n\n");
  printf("testing fn_flux.cpp function\n\n");

  std::vector<double> fn = fn_flux(ul, beta);
  PrintingContainer(fn);
  
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
  
  

  printf("#######################################################################\n\n");

  return 0;
}
