#include <math.h>
#include "diffusion_tensor.h"
#include "diffusion_matrix.h"
#include "utilities.h"

void diffusion_matrix(std::vector<std::vector<double>>& u, double h, double a, std::map<int, std::vector<std::vector<double>>> &D, std::map<int, std::vector<std::vector<double>>> &L, std::map<int, std::vector<std::vector<double>>> &U)
{
	int m = u.size(); int n = u[0].size();

	double f = -a * 0.5 / (h * h);
	std::map<int, std::vector<std::vector<double>>> B = diffusion_tensor(u);
	std::vector<std::vector<double>> Mele(m , std::vector<double>(m));

	for (int i = 0; i < n; i++) {
		D[i] = Mele; L[i] = Mele; U[i] = Mele;
	}

	double z;
	for (int p = 0; p < m; p++) 
	{
		for (int j = 0; j < n; j++) 
		{
			if (j > 0)
			{
				for (int i = 0; i < m; i++) 
				{
					//K(p, j)=K(p, j)-(B(p, i, j)+B(p, i, j-1))*(u(i, j)-u(i, j-1));
					z = B[j][p][i] + B[j - 1][p][i];
					D[j][p][i] = D[j][p][i] - f * z;
					L[j - 1][p][i] = L[j - 1][p][i] + f * z;
				}
			}
			if (j < n - 1)
			{
				for (int i = 0; i < m; i++) 
				{
					//K(p, j)=K(p, j)+(B(p, i, j+1)+B(p, i, j))*(u(i, j+1)-u(i, j));
					z = B[j + 1][p][i] + B[j][p][i];
					U[j][p][i] = U[j][p][i] + f * z;
					D[j][p][i] = D[j][p][i] - f * z;
				}
			}
		}
	}

	for (int j = 0; j < n; j++) 
		for (int i = 0; i < m; i++)
			D[j][i][i] = D[j][i][i] + 1;
}