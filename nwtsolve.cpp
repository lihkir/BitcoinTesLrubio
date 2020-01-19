#include <stdio.h>
#include <math.h>
#include "bc.h"
#include "diffusion_matrix.h"
#include "newton_matrix.h"
#include "residual.h"
#include "trdsolve.h"
#include "nwtsolve.h"
#include "test_cases.h"
#include "update_utilities.h"

void nwtsolve(double a, Matrix<double> &b, Matrix<double> &u, double h, double dt, int maxits, double tol)
{
    struct test_cases* pt_test = get_tests();

    int m = size(u, 1);
    int n = size(u, 2) - 2*pt_test->gc;

    Matrix<double> *ptv = new Matrix<double>(m, n); Matrix<double> &v = *ptv;
    update_inside(v, u, pt_test->gc);

    Block<double> *ptD = new Block<double>(n, m, m); Block<double> &D = *ptD;
    Block<double> *ptL = new Block<double>(n, m, m); Block<double> &L = *ptL;
    Block<double> *ptU = new Block<double>(n, m, m); Block<double> &U = *ptU;

    Matrix<double> *ptr = new Matrix<double>(m, n); Matrix<double> &r = *ptr;
    Matrix<double> *ptrh = new Matrix<double>(m, n); Matrix<double> &rh = *ptrh;
    Matrix<double> *ptdu = new Matrix<double>(m, n); Matrix<double> &du = *ptdu;
    
    diffusion_matrix(v, h, a, D, L, U);
    update_inside(r, b, 0);
    residual(D, L, U, v, r);

    Matrix<double> *ptbj = new Matrix<double>(m, 1); Matrix<double> &bj = *ptbj;
    Matrix<double> *ptrhj = new Matrix<double>(m, 1); Matrix<double> &rhj = *ptrhj;
    Matrix<double> *ptrj = new Matrix<double>(m, 1); Matrix<double> &rj = *ptrj;

    /** r=b-(L|D|U)*v **/    
    double nrm2_b = 0;
    for (int j = 1; j <= n; j++)
    {
        get_col(bj, b, j);
        nrm2_b += dot_product(bj, bj);
    }

    double alpha, phi; 

    for (int its = 1; its <= maxits; its++)
    {
        newton_matrix(v, h, a, D, L, U);

        update_inside(du, r, 0);
        trdsolve(D, L, U, du);

        update_inside(rh, r, 0);
        residual(D, L, U, du, rh);
  
        double nrm2_r0 = 0;
        for (int j = 1; j <= n; j++)
        {
            get_col(rhj, rh, j);
            nrm2_r0 += dot_product(rhj, rhj);
        }
            
        diffusion_matrix(v, h, a, D, L, U);

        update_inside(r, b, 0);
        residual(D, L, U, v, r);
  
        double phi0 = 0;
        for (int j = 1; j <= n; j++)
        {   
            get_col(rj, r, j);
            phi0 += dot_product(rj, rj);
        }
        
        alpha = 1;
        phi = 2*phi0 + 1;

        while (alpha > 1e-4 && phi>=phi0)
        {
            for (int j = 1; j <= n; j++)
                for (int i = 1; i <= m; i++)
                    v(i, j) = v(i, j) + alpha*du(i, j);
            
            diffusion_matrix(v, h, a, D, L, U);

            update_inside(r, b, 0);
            residual(D, L, U, v, r);
    
            phi = 0;
            for (int j = 1; j <= n; j++)
            {
                get_col(rj, r, j);
                phi += dot_product(rj, rj);
            }

            alpha = alpha/2;
        }
        /** calcule r=b-B(v)*v **/
        diffusion_matrix(v, h, a, D, L, U);

        update_inside(r, b, 0);
        residual(D, L, U, v, r);

        double nrm2_r=0;
        for (int j = 1; j <= n; j++)
        {
            get_col(rj, r, j);
            nrm2_r += dot_product(rj, rj);   
        }
    
        if (nrm2_r > 0)
            printf("\t%3d %e %e\n", its, sqrt(nrm2_r), sqrt(nrm2_r0));

        if (nrm2_r < tol*tol*nrm2_b)
        {
            if ( nrm2_r > 0 )
            {
                printf("\n----------\n");
                break;
            }
        }
    }

    update_inside(u, v, pt_test->gc);
    bc(u);

    delete ptbj;
    delete ptD;
    delete ptdu;
    delete ptL;
    delete ptr;
    delete ptrh;
    delete ptrhj;
    delete ptrj;
    delete ptU;
    delete ptv;
}