#include <stdio.h>
#include <math.h>
#include <assert.h> 
#include <map>
#include "diffusion_tensor.h"
#include "test_cases.h"
#include "utilities.h"
#include "sed_hsf.h"
#include "sed_hsf_der.h"
#include "solid_stress.h"
#include "solid_stress_der.h"
#include "solid_stress_der_der.h"
#include "deltak.h"

std::map<int, std::vector<std::vector<double>>> der_diffusion_tensor(std::vector<std::vector<double>>& u, int r)
{
  struct test_cases* pt_test = get_tests();
  
  int m = u.size(); int n = u[0].size();
  
  std::map<int, std::vector<std::vector<double>>> Bsol;
  std::vector<std::vector<double>> Bele(m, std::vector<double>(n));
  
  double dphi, phit, wphi, dwphi, wopq, sige, sigd, sigdd, pwpsq, cste_one, wophi, dp1mp, p1mpq, phidelta, cste_two, cste_new;

  for (int k = 0; k < n; k++)
  {
    std::vector<double> phi = col(u, k); 
    
    dphi  = dot_product(pt_test->delta, phi);
    phit  = vector_sum(phi);
    wphi  = sed_hsf(phit);
    dwphi = sed_hsf_der(phit);
    sige  = solid_stress(phit);
    sigd  = solid_stress_der(phit);
    sigdd = solid_stress_der_der(phit);
    pwpsq = (phit*dwphi - wphi)/pow(phit, 2);
    wophi = wphi/phit; 
    dp1mp = (1 - 2*phit)/pow(phit*(1 - phit), 2); 
    p1mpq = 1/phit*(1 - phit);
    
    if (phit > 1) throw std::invalid_argument("\nError: phit > 1 !!\n"); 
    Bsol[k] = Bele;
    
    for (int i = 0; i < m; i++)
	  {
      phidelta = phi[i]*(pt_test->delta[i] - dphi);
	    for (int j = 0; j < m; j++)
	    {
        cste_new = pt_test->delta[i]*deltak(i, j) - pt_test->delta[j]*phi[i] - phidelta*p1mpq;
        cste_one = phidelta*sigd - cste_new*sige;
        cste_two = phidelta*sigdd - (dp1mp*phidelta*sige + sigd*cste_new);
	      Bsol[k][i][j] = pt_test->muog*(pwpsq*cste_one + wophi*cste_two);
	    }
	  }
  }
  return Bsol;
}