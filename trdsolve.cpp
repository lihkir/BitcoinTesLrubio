#include <stdio.h>
#include <math.h>
#include "trdsolve.h"
#include "utilities.h"
#include "update_utilities.h"

void trdsolve(std::map<int, std::vector<std::vector<double>>>& A1, std::map<int, std::vector<std::vector<double>>>& B1, std::map<int, std::vector<std::vector<double>>>& C1, std::vector<std::vector<double>>& b)
{
    int n = A1.size();
    int m = A1[0].size();

    std::map<int, std::vector<std::vector<double>>> A = A1;
    std::map<int, std::vector<std::vector<double>>> B = B1;
    std::map<int, std::vector<std::vector<double>>> C = C1;

    lutrb(A, B, C);

    std::vector<std::vector<double>> bz = sub_col(b, 0);
    fwdsolve(A[0], bz); UpdateCol(b, bz,0);
    
    for (int i = 1; i < n; i++) 
    {
        std::vector<std::vector<double>> bl = sub_col(b, i - 1);
        std::vector<std::vector<double>> bc = sub_col(b, i);
        maxpy(B[i-1], bl, bc);  
        fwdsolve(A[i], bc); UpdateCol(b, bc, i);
    }

    std::vector<std::vector<double>> bn = sub_col(b, n-1);
    bwdsolve(A[n-1], bn); UpdateCol(b, bn, n-1);

    for (int i = n - 2; i >= (0); i -= 1) 
    {
        std::vector<std::vector<double>> br = sub_col(b, i + 1);
        std::vector<std::vector<double>> bc = sub_col(b, i);
        maxpy(C[i], br, bc); 
        bwdsolve(A[i], bc); UpdateCol(b, bc, i);
    }
}

void lutrb(std::map<int, std::vector<std::vector<double>>> &A, std::map<int, std::vector<std::vector<double>>> &B, std::map<int, std::vector<std::vector<double>>> &C)
{   
    int n = A.size();    
    for (int i = 0; i < n-1; i++) 
    {
        nplu(A[i]);
        utrsolve(A[i], B[i]);
        fwdsolve(A[i], C[i]);
        maxpy(B[i], C[i], A[i+1]);
    }
    nplu(A[n-1]);
}

void nplu(std::vector<std::vector<double>> &A)
{  
    int n = A.size();
    for (int k = 0; k < n; k++)
    {
        if (fabs(A[k][k]) < 1e-15) printf("\nPivot too small\n");
        for (int i = k+1; i < n; i++)
        {
            A[i][k] = A[i][k]/A[k][k];
            for (int j = k+1; j < n; j++)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
        }
    } 
}

void fwdsolve(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
    int s = B[0].size();
    int n = A.size();
    for (int k = 0; k < s; k++) 
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < i-1; j++) 
                B[i][k] = B[i][k] - A[i][j] * B[j][k];
}

void bwdsolve(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
    int s = B[0].size();
    int n = A.size();
    for (int k = 0; k < s; k++) 
    {
        for (int i = n - 1; i >= (0); i-=1) 
        {
            for (int j = i + 1; j < n; j++) B[i][k] = B[i][k] - A[i][j] * B[j][k];
            B[i][k] = B[i][k]/A[i][i];
        }
    }
}

void utrsolve(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
    int s = B[0].size();
    int n = A.size();
    for (int k = 0; k < s; k++) 
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < i-1; j++) B[k][i] = B[k][i] - A[j][i] * B[k][j];
            B[k][i] = B[k][i]/A[i][i];
        }
    }
}

void maxpy(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &y)
{
    int m = y.size();
    int l = y[0].size();
    int n = A[0].size();
 
    for (int k = 0; k < l; k++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i][k] = y[i][k] - A[i][j] * x[j][k];
}