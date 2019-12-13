#include <stdio.h>
#include <math.h>
#include <assert.h> 
#include <map>
#include "test_cases.h"
#include "utilities.h"
#include "diffusion_tensor.h"
#include "sed_hsf.h"
#include "solid_stress.h"
#include "solid_stress_der.h"
#include "deltak.h"

std::map<int, std::vector<std::vector<double>>> diffusion_tensor(std::vector<std::vector<double>>& u)
{
  struct test_cases* pt_test = get_tests();
  
  int m = u.size();
  int n = u[0].size();
  
  std::map<int, std::vector<std::vector<double>>> Bsol;
  std::vector<std::vector<double>> Bele(m, std::vector<double>(m));
  
  double dphi, phit, wphi, wovq, sige, sigd, q1mp, phidelta, cste_new;

  for (int k = 0; k < n; k++)
  {	
    std::vector<double> phi = col(u, k); 
    
    dphi = dot_product(pt_test->delta, phi);
    phit = vector_sum(phi);
    wphi = sed_hsf(phit);
    wovq = wphi/phit;
    q1mp = 1/phit*(1 - phit);
    sige = solid_stress(phit);
    sigd = solid_stress_der(phit);

    if (phit > 1) throw std::invalid_argument("\nError: phit > 1 !!\n");
    Bsol[k] = Bele;
        
    for (int i = 0; i < m; i++)
    {
      phidelta = phi[i]*(pt_test->delta[i] - dphi);
      for (int j = 0; j < m; j++)
      {
        cste_new = pt_test->delta[i]*deltak(i, j) - pt_test->delta[j]*phi[i] - q1mp*phidelta;
        Bsol[k][i][j] = pt_test->muog*wovq*(phidelta*sigd - cste_new*sige); 
      }
    }
  }
  return Bsol;
}