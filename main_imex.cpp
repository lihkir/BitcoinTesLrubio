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
  	double T = atof(argv[8]);
	global::diff_idx = atoi(argv[9]);
	const char *outname= argv[10];
	
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

	clock_t start, end;
  	start = clock();
	rkimex(A, b, Ah, bh, u0, T, cfl, pt_test->L, pt_test->convec_type, imex_type);
	end = clock();
  	double CPUTIME = double(end - start) / CLOCKS_PER_SEC;
	
  	printf("\nCPUTIME rkimex T=%.15f=> %.15f\n", T, CPUTIME);

	ofstream fileout(outname);
  	for (int i = 0; i < pt_test->M_rows; i++)
    	for (int j = 0; j < N_cols; j++) 
      		fileout << std::fixed << std::setprecision(16) << std::scientific << u0[i][j] << endl;
  	fileout.close();

  	delete pu0;
	delete pt_test->delta;
	delete ptA;
	delete ptb;
	delete ptAh;
	delete ptbh;

	return 0;
}