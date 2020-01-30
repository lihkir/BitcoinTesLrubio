#include <fstream>
#include <cmath>
#include <string.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <assert.h>
#include "test_cases.h"
#include "containers.h"
#include "rkimex.h"

using namespace std;

namespace global
{ 
  extern int idx_test;
  extern int idx_q;
  extern double D0;
  extern int int_form;
  extern int diff_idx;
}

Vector<double> parse_times(char * str_times)
{
  char * p;
  int commas = 0;
  for (p = str_times; *p; p++) 
    if (*p == ',') commas++;
  
  Vector<double> Ta(commas + 1);
  
  Ta(1) = atof(strtok(str_times, ","));
  int count = 0;
  while ((p = strtok(NULL, ",")))
    Ta(++count + 1) = atof(p);
  
  return Ta;  
}

int main(int argc, char* argv[])
{
  printf("\nStarting pvm_sed\n\n");
  /** M_rows := number of species and N_cols := number of nodes **/
  
  int N_cols = atoi(argv[1]);	
  global::idx_test = atoi(argv[2]);
  global::D0 = atof(argv[3]);
  global::idx_q = atoi(argv[4]);
  global::int_form = atoi(argv[5]);
  double cfl = atof(argv[6]);
  int imex_type = atoi(argv[7]);
  Vector<double> Ta = parse_times(argv[8]);
  global::diff_idx = atoi(argv[9]);
  
  struct test_cases* pt_test = get_tests();	
  
  Matrix<double> *pu0 = new Matrix<double>(pt_test->M_rows, N_cols); Matrix<double> &u0 = *pu0;
  Matrix<double> &delta = *pt_test->delta;
  
#include "init.h"
  
  Matrix<double> *ptA = new Matrix<double>(1, 1); Matrix<double> &A = *ptA;
  Vector<double> *ptb = new Vector<double>(1); Vector<double> &b = *ptb;
  A(1, 1) = 0.5;
  b(1) = 1;
  
  Matrix<double> *ptAh = new Matrix<double>(2, 2); Matrix<double> &Ah = *ptAh;
  Vector<double> *ptbh = new Vector<double>(2); Vector<double> &bh = *ptbh;
  Ah(2, 1) = 1;
  bh(1) = 0.5;
  bh(2) = 0.5;
  
  print2D(Ah, 8, "Ah");
  print1D(bh, 8, "bh");
  print2D(A, 8, "A");
  print1D(b, 8, "b");
  
  printf("\n##########################################\n\n");
  rkimex(A, b, Ah, bh, u0, Ta, cfl, pt_test->L, pt_test->convec_type, imex_type);
  
  delete pu0;
  delete pt_test->delta;
  delete ptA;
  delete ptb;
  delete ptAh;
  delete ptbh;
  
  return 0;
}