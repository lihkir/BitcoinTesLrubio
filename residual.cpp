#include <stdio.h>
#include <math.h>
#include "residual.h"

void residual(std::map<int, std::vector<std::vector<double>>>& D, std::map<int, std::vector<std::vector<double>>>& L, std::map<int, std::vector<std::vector<double>>>& U, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& r)
{  
    int m = r.size();
    int n = D.size();
    
    for (int i = 0; i < m; i++) 
    {
        for (int k = 0; k < m; k++)
            r[i][0] = r[i][0] - D[0][i][k] * v[k][0] - U[0][i][k] * v[k][1];

        for (int j = 1; j < n-1; j++ ) 
            for (int k = 0; k < m; k++) 
                r[i][j] = r[i][j] - L[j-1][i][k] * v[k][j-1] - D[j][i][k] * v[k][j] - U[j][i][k] * v[k][j+1];                   

        for (int k = 0; k < m; k++)
            r[i][n-1] = r[i][n-1] - L[n-2][i][k] * v[k][n-2] - D[n-1][i][k] * v[k][n-1];
    }
}