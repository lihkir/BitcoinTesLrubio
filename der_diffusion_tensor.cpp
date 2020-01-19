#include "diffusion_tensor.h"
#include "containers.h"
#include "test_cases.h"
#include "solid_stress_der_der.h"
#include "solid_stress_der.h"
#include "solid_stress.h"
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "kronecker.h"
#include "update_utilities.h"

namespace global 
{
	extern double D0; 
    extern int diff_idx;
}

void der_diffusion_tensor(Matrix<double>& u, Block<double> &B, int r)
{
	struct test_cases* pt_test = get_tests();

	int m = size(u, 1); 
	int n = size(u, 2);

	Matrix<double> &delta = *pt_test->delta;
	Matrix<double> *ptphik =  new Matrix<double>(m, 1); Matrix<double> &phik = *ptphik;
	
	double phit, p2, quot, sigedd, siged, sige, wphit, wphitd, theta, psih;

	for (int k = 1; k <= n; k++)
	{
		get_col(phik, u, k);
		phit = sum1D(phik);
		p2 = dot_product(delta, phik);

		if (phit > 1)
    	{
      		printf("\nFor column %d\n", k);
      		throw std::invalid_argument("\nError inside der_diffusion_tensor: phit > 1 !!\n");
    	}

    	if (global::diff_idx != 0)
		{
			quot = 1/phit*(1 - phit);
			sigedd = solid_stress_der_der(phit);
        	siged  = solid_stress_der(phit);
        	sige   = solid_stress(phit);
        	wphit  = sed_hsf(phit);
        	wphitd = sed_hsf_der(phit);

			for (int i = 1; i <= m; i++)
			{
				theta = phik(i)*(delta(i) - p2);
				for (int j = 1; j <= m; j++)
				{
					psih = delta(i)*kronecker(i, j) - delta(j)*phik(i) - theta/phit;
                	B(k, i, j) = pt_test->mu_g*((phit*wphitd - wphit)*(theta*siged - psih*sige*quot)/pow(phit, 2) 
							   + wphit*(theta*sigedd - (theta*sige/pow(phit, 2) + psih*((1 - phit)*siged + sige)/(1 - phit))/(1 - phit))/phit);					
				}			
			}
		}
		else
		{
			for (int i = 1; i <= m; i++) 
				B(k, i, i) = -global::D0*pt_test->nexp*pow(1 - phit, pt_test->nexp - 1);
		}
	}
	
	delete ptphik;
	delete pt_test->delta;
	delete pt_test;
}