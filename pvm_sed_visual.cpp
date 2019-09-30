#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <assert.h>
#include "test_cases.h"
#include "utilities.h"
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "average.h"
#include "bc.h"
#include "fn_flux.h"
#include "matmult.h"
#include "hornerm.h"
#include "lminmax_charspeed.h"
#include "minmax_charspeed.h"
#include "diffusion_tensor.h"
#include "der_diffusion_tensor.h"
#include "charspeed.h"
#include "apply_diffus.h"
#include "diffus.h"
#include "diffus_charspeed.h"

using namespace std;

namespace global
{ 
	extern int idx_test;
	extern double D0;
}

int main(int argc, char* argv[])
{
	printf("\n\nStarting pvm_sed\n\n");
	/** M_rows := number of species and N_cols := number of nodes **/

	int N_cols = atoi(argv[1]);	
	global::idx_test = atoi(argv[2]);
	global::D0 = atof(argv[3]);

	srand(unsigned int(time(NULL)));
	int Prec = 16;

	printf("\n#######################################################################\n");
	printf("\ntesting: beta_cases.cpp function\n\n");

	struct test_cases* pt_test = get_tests();
	printf("D0 = %0.16f\t pt_test->nexp = %f\n\n", global::D0 , pt_test->nexp);
	PrintingContainer(pt_test->delta, Prec);

	printf("#######################################################################\n\n");

	double sed_hsfd = sed_hsf(0.5);
	printf("testing: sed_hsfd = %0.16f\n\n", sed_hsfd);

	double sed_hsfd_der = sed_hsf_der(1);
	printf("testing: sed_hsfd = %0.16f\n\n", sed_hsfd_der);

	printf("\n#######################################################################\n");
	printf("\ntesting: average.cpp function\n\n");

	std::vector<double> ul = RandomVector<double>(pt_test->M_rows, 1e-14);
	std::vector<double> ur = RandomVector<double>(pt_test->M_rows, 1e-14);
	std::vector<double> az = RandomVector<double>(pt_test->M_rows, 1e-14);

	PrintingContainer(ul, Prec);
	PrintingContainer(ur, Prec);
	std::vector<double> ua = average(ul, ur, 0.5);
	PrintingContainer(ua, Prec);

	printf("\n#######################################################################\n");
	printf("\ntesting: bc.cpp function\n\n");

	std::vector<std::vector<double>> m   = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> z   = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> p   = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> A	 = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> x   = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> pAx = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> pr  = RandomMatrix<double>(pt_test->M_rows, N_cols, 1e-14);
	std::vector<std::vector<double>> Sl  = RandomMatrix<double>(2, N_cols, 1e-14);

	printf("\nm, z, p matrices before ghost:\n\n");
	PrintingContainer(m, Prec);
	PrintingContainer(z, Prec);
	PrintingContainer(p, Prec);
	bc(m);
	bc(z);
	bc(p);
	bc(Sl);
	printf("\nm, z, p matrices after ghost:\n\n");
	PrintingContainer(m, Prec);
	PrintingContainer(z, Prec);
	PrintingContainer(p, Prec);

	printf("\n#######################################################################\n");

	printf("testing fn_flux.cpp function\n\n");

	std::vector<double> fn = fn_flux(ul);
	int Prec32 = 32;
	PrintingContainer(fn, Prec32);

	printf("#######################################################################\n\n");
	printf("\ntesting hornerm.cpp\n");

	bc(A);
	bc(x);
	bc(pr);
	bc(pAx);
	printf("\naz, A, x, pr, pAx: before\n");
	PrintingContainer(az, Prec);
	PrintingContainer(A, Prec);
	PrintingContainer(x, Prec);
	PrintingContainer(pr, Prec);
	PrintingContainer(pAx, Prec);
	hornerm(A, x, az, pr, pAx);
	printf("\naz, A, x, pr, pAx: later\n");
	PrintingContainer(az, Prec);
	PrintingContainer(A, Prec);
	PrintingContainer(x, Prec);
	PrintingContainer(pr, Prec32);
	PrintingContainer(pAx, Prec32);

	printf("#######################################################################\n\n");
	printf("testing matmult.cpp function\n\n");

	printf("\nm, z, p before:\n");
	PrintingContainer(m, Prec);
	PrintingContainer(z, Prec);
	PrintingContainer(p, Prec);
	matmult(m, z, p);
	printf("\nm, z, p later:\n");
	PrintingContainer(m, Prec);
	PrintingContainer(z, Prec);
	PrintingContainer(p, Prec);

	printf("#######################################################################\n\n");
	printf("testing lminmax_charspeed.cpp function\n\n");

	std::vector<double> Sleft = lminmax_charspeed(ul, ur);
	PrintingContainer(Sleft, Prec);

	printf("#######################################################################\n\n");
	printf("testing minmax_charspeed.cpp function\n\n");

	std::vector<std::vector<double>> Sminmax = minmax_charspeed(m);
	PrintingContainer(Sminmax, Prec);

	printf("#######################################################################\n\n");
	printf("testing diffusion_tensor.cpp function\n\n");

	std::map<int, std::vector<std::vector<double>>> MapDiff = diffusion_tensor(m);
	for (unsigned int i = 0; i < m.size(); i++) { PrintingContainer(MapDiff[i], Prec); }

	printf("#######################################################################\n\n");
	printf("testing der_diffusion_tensor.cpp function\n\n");

	std::map<int, std::vector<std::vector<double>>> MapDerDiff = der_diffusion_tensor(m);
	for (unsigned int i = 0; i < m.size(); i++) { PrintingContainer(MapDerDiff[i], Prec); }

	printf("#######################################################################\n\n");
	printf("testing charspeed.cpp function\n\n");

	double chrs = charspeed(m);
	printf("charspeed = %0.16f\n\n", chrs);

	printf("#######################################################################\n\n");
	printf("testing apply_diffus.cpp function\n\n");

	std::vector<std::vector<double>> K = apply_diffus(m, z, 0.5);
	PrintingContainer(K, 32);

	printf("#######################################################################\n\n");
	printf("testing diffus.cpp function\n\n");

	std::vector<std::vector<double>> Kd = diffus(m,0.5);
	PrintingContainer(Kd, 32);

	printf("#######################################################################\n\n");
	printf("testing diffus_charspeed.cpp function\n\n");

	printf("diffus_charspeed(m) = %0.16f\n\n", diffus_charspeed(m));

	printf("#######################################################################\n\n");

	std::cin.get();
	return 0;
}