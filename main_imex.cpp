#include <fstream>
#include <cmath>
#include <string.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <assert.h>
#include <iomanip>
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
#include "diffusion_matrix.h"
#include "jacobiana.h"
#include "jacobiana_dec.h"
#include "weno5.h"
#include "convec.h"
#include "residual.h"
#include "lsolve.h"
#include "rkimex.h"
#include "solid_stress.h"
#include "solid_stress_der.h"
#include "deltak.h"

using namespace std;

namespace global
{ 
  extern int idx_test;
  extern int idx_q;
  extern double D0;
  extern int int_form;
  extern double beta;
}

std::vector<double> parse_times(char * str_times)
{
	char * p;
	int commas=0;
	for ( p=str_times; *p; p++ ) 
		if ( *p == ',' ) commas++;

	std::vector<double> Ta(commas+1);

  	Ta[0] = atof(strtok(str_times, ","));
	int count = 0;
	while ((p = strtok(NULL, ",")))
		Ta[++count] = atof(p);
	
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
	std::vector<double> Ta = parse_times(argv[8]);
	const char *outname= argv[9];

 	struct test_cases* pt_test = get_tests();
  
  	std::vector<std::vector<double>> u0(pt_test->M_rows, std::vector<double>(N_cols));
  	global::beta = 1e-15;

#include "init.h"

	std::vector<std::vector<double>> A;
  	std::vector<double> b;
  	std::vector<std::vector<double>> Ah;
  	std::vector<double> bh;
  
	if (imex_type == 0) 
    { /** NI **/
		std::vector<double> Arow; Arow.push_back(0.5);
      	A.push_back(Arow);
      	b.push_back(1);
      
      	std::vector<double> Ahrow; Ahrow.push_back(0); Ahrow.push_back(0.5);
      	Ah.push_back(std::vector<double>(2)); Ah.push_back(Ahrow);
      	bh.push_back(0); bh.push_back(1);
    }
	else
    {
    	std::vector<double> Arow; Arow.push_back(0.5); Arow.push_back(0.5);
      	A.push_back(std::vector<double>(2)); A.push_back(Arow);
      	b.push_back(0.5); b.push_back(0.5);
      
      	std::vector<double> Ahrow; Ahrow.push_back(1); Ahrow.push_back(0);
      	Ah.push_back(std::vector<double>(2)); Ah.push_back(Ahrow);
      	bh.push_back(0.5); bh.push_back(0.5);
    }
  	printf("\nFunction: rkimex.cpp\n");
	rkimex(A, b, Ah, bh, u0, Ta, cfl, pt_test->L, pt_test->convec_type, imex_type, global::idx_test);
  	
	ofstream fileout(outname);
  	for (int i = 0; i < pt_test->M_rows; i++)
    	for (int j = 0; j < N_cols; j++) 
      		fileout << std::fixed << std::setprecision(16) << std::scientific << u0[i][j] << endl;
  	fileout.close();

	/****************************************************************************************/
	printf("\nTesting solid_stress function:\n");
	printf("\nsolid_stress(phi) = %0.16f\n", solid_stress(0.5));

	printf("\nTesting solid_stress function derivative:\n");
	printf("\nsolid_stress_der(phi) = %0.16f\n", solid_stress_der(0.5));

	printf("\nTesting deltak function:\n");
	printf("\ndeltak(i, j) = %0.16f\n", deltak(2, 2));
	
	printf("\nTesting sed_hsf function:\n");
	printf("\nsed_hsf(d) = %0.16f\n", sed_hsf(0.4539393939393940));

	printf("\nTesting sed_hsf_der function:\n");
	printf("\nsed_hsf_der(d) = %0.16f\n", sed_hsf_der(0.4539393939393940));
	
	/****************************************************************************************/
	return 0;
}