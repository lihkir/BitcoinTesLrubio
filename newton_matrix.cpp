#include <stdio.h>
#include <math.h>
#include "der_diffusion_tensor.h"
#include "diffusion_tensor.h"
#include "newton_matrix.h"

void newton_matrix(std::vector<std::vector<double>> &u, double h, double a, std::map<int, std::vector<std::vector<double>>> &D, std::map<int, std::vector<std::vector<double>>> &L, std::map<int, std::vector<std::vector<double>>> &U)
{
    int n = u[0].size();
    int m = u.size();

    double z;
    double f = -a*0.5/(h*h);
    std::map<int, std::vector<std::vector<double>>> B = diffusion_tensor(u);

	std::vector<std::vector<double>> Mele(m , std::vector<double>(m));
	for (int i = 0; i < n; i++) {
		D[i] = Mele; L[i] = Mele; U[i] = Mele;
	}

    for (int r = 0; r < m; r++)
    {
        std::map<int, std::vector<std::vector<double>>> Bp = der_diffusion_tensor(u, r);
        for (int p = 0; p < m; p++) 
        {
            for (int i = 1; i < n; i++) 
            {
                z = B[i][p][r] + B[i-1][p][r];
                D[i][p][r] = D[i][p][r] - f*z;
                L[i-1][p][r] = L[i-1][p][r] + f*z;
                
                for (int q = 0; q < m; q++)
                    D[i][p][r] = D[i][p][r] - f*Bp[i][p][q]*(u[q][i] - u[q][i-1]);

                for (int q = 0; q < m; q++)
                    L[i-1][p][r] = L[i-1][p][r] - f*Bp[i-1][p][q]*(u[q][i] - u[q][i-1]);
            }
            for (int i = 0; i < n-1; i++) 
            {
                z = B[i+1][p][r] + B[i][p][r];
                U[i][p][r] = U[i][p][r] + f*z;
                D[i][p][r] = D[i][p][r] - f*z;

                for (int q = 0; q < m; q++) 
                    D[i][p][r] = D[i][p][r] + f*Bp[i][p][r]*(u[q][i+1] - u[q][i]);
                
                for (int q = 0; q < m; q++)
                    U[i][p][r] = U[i][p][r] + f*Bp[i+1][p][q]*(u[q][i+1] - u[q][i]);
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int r = 0; r < m; r++)
            D[i][r][r] = D[i][r][r] + 1;
}