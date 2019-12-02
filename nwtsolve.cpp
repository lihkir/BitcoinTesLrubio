#include <stdio.h>
#include <math.h>
#include "bc.h"
#include "diffusion_matrix.h"
#include "newton_matrix.h"
#include "residual.h"
#include "trdsolve.h"
#include "nwtsolve.h"
#include "test_cases.h"
#include "utilities.h"
#include "update_utilities.h"

void nwtsolve(double a, std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &u, double h, double dt, int maxits, double tol)
{
    struct test_cases* pt_test = get_tests();

    int m = u.size();
    int n = u[0].size() - 2*pt_test->gc;

    std::vector<std::vector<double>> v = sub_matrix(u, pt_test->gc);
    
    std::map<int, std::vector<std::vector<double>>> D;
    std::map<int, std::vector<std::vector<double>>> L;
    std::map<int, std::vector<std::vector<double>>> U;

    diffusion_matrix(v, h, a, D, L, U);
    std::vector<std::vector<double>> r = sub_matrix(b, 0);
    residual(D, L, U, v, r); /** r = b-(L||D||U)*v **/
    
    double nrm2_b = 0;
    for (int j = 0; j < n; j++) nrm2_b += square_sum(sub_vector(b, j));

    for (int its = 0; its < maxits; its++) 
    {
        newton_matrix(v, h, a, D, L, U);
        std::vector<std::vector<double>> du = sub_matrix(r, 0);
        trdsolve(D, L, U, du);
        std::vector<std::vector<double>> r1 = sub_matrix(r, 0);
        residual(D, L, U, du, r1);
        
        double nrm2_r0 = 0;
        for (int j = 0; j < n; j++) nrm2_r0 += square_sum(sub_vector(r1, j));

        diffusion_matrix(v, h, a, D, L, U);

        std::vector<std::vector<double>> r = sub_matrix(b, 0);
        residual(D, L, U, v, r);
        
        double phi0 = 0;
        for (int j = 0; j < n; j++) phi0 =+ square_sum(sub_vector(r, j));
        
        double alpha = 1;
        double phi = 2*phi0 + 1;

        while (alpha > 1e-4 && phi >= phi0) 
        {
            for (int j = 0; j < n; j++) 
                for (int i = 0; i < m; i++) 
                    v[i][j] = v[i][j] + alpha*du[i][j];
            
            diffusion_matrix(v, h, a, D, L, U);
            std::vector<std::vector<double>> r = sub_matrix(b, 0);
            residual(D, L, U, v, r);
            
            phi = 0;
            for (int j = 0; j < n; j++) phi += square_sum(sub_vector(r, j)); 
            alpha = alpha/2;
        }
  
        diffusion_matrix(v, h, a, D, L, U);
        std::vector<std::vector<double>> rh = sub_matrix(b, 0);
        residual(D, L, U, v, rh);
        
        double nrm2_r=0;
        for (int j = 0; j < n; j++) nrm2_r += square_sum(sub_vector(rh, j)); 
        
        double beta = nrm2_r;
        if  (nrm2_r > 0) printf("\t%3d %e %e\n", its, sqrt(nrm2_r), sqrt(nrm2_r0));
  
        if  (nrm2_r < tol*tol*nrm2_b) 
        {
            if  ( nrm2_r > 0 ) printf("----------\n");
            break;
        }
    }

    for (unsigned int i = 0; i < v.size(); i++)
        for (unsigned int j = 0; j < v[0].size(); j++)
            u[i][pt_test->gc + j] = v[i][j];
    bc(u);
}