#include "trdsolve.h"
#include "update_utilities.h"

void trdsolve(Block<double> &A1, Block<double> &B1, Block<double> &C1, Matrix<double> &b)
{
    int n = size(A1, 1);
    int m = size(A1, 2);

    Block<double> A(n, m, m);
    Block<double> B(n, m, m);
    Block<double> C(n, m, m);
    
    copy_block(A, A1);
    copy_block(B, B1);
    copy_block(C, C1);
    
    lutrb(A, B, C);

    Matrix<double> bf(m, 1);
    Matrix<double> bl(m, 1);
    Matrix<double> br(m, 1);
    Matrix<double> bc(m, 1);
    Matrix<double> be(m, 1);

    get_col(bf, b, 1);
    fwdsolve(A[0], bf);
    update_col(b, bf, 1);
    
    for (int i = 2; i <= n; i++)
    {
        get_col(bl, b, i - 1);
        get_col(bc, b, i);
        
        maxpy(B[i - 2], bl, bc);
        fwdsolve(A[i - 1], bc);

        update_col(b, bl, i - 1);
        update_col(b, bc, i);
    }

    get_col(be, b, n);
    bwdsolve(A[n - 1], be);
    update_col(b, be, n);

    for (int i = n - 1; i >= (1); i -= 1)
    {
        get_col(br, b, i + 1);
        get_col(bc, b, i);

        maxpy(C[i - 1], br, bc);
        bwdsolve(A[i - 1], bc);

        update_col(b, br, i + 1);
        update_col(b, bc, i);
    }
}

void lutrb(Block<double> &A, Block<double> &B, Block<double> &C)
{
    int n = size(A, 1);

    for (int i = 1; i <= n-1; i++)
    {
        nplu(A[i - 1]);
        utrsolve(A[i - 1], B[i - 1]);
        fwdsolve(A[i - 1], C[i - 1]);
        maxpy(B[i - 1], C[i - 1], A[i]);
    }
    nplu(A[n - 1]);
}

void nplu(Matrix<double> &A)
{
    int n = size(A, 1);
    
    for (int k = 1; k <= n; k++)
    {
        if (abs(A(k,k)) < 1e-15)
            throw std::invalid_argument("\nPivot too small!!\n");
        
        for (int i = k + 1; i <= n; i++)
        {
            A(i, k) = A(i, k)/A(k, k);
            for (int j = k + 1; j <= n; j++)
                A(i, j) = A(i, j) - A(i, k)*A(k, j);
        }
    }
}

void fwdsolve(Matrix<double> &A, Matrix<double> &B)
{
    int s = size(B, 2);
    int n = size(A, 1);

    for (int k = 1; k <= s; k++)
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= i - 1; j++)
                B(i, k) = B(i, k) - A(i, j)*B(j, k);
}

void bwdsolve(Matrix<double> &A, Matrix<double> &B)
{
    int s = size(B, 2);
    int n = size(A, 1);
    
    for (int k = 1; k <= s; k++)
    {
        for (int i = n; i >= (1); i -= 1)
        {
            for (int j = i + 1; j <= n; j++)
                B(i, k) = B(i, k) - A(i, j)*B(j, k);
            B(i, k) = B(i, k)/A(i, i);
        }
    }
}

void utrsolve(Matrix<double> &A, Matrix<double> &B)
{
    int s = size(B, 2);
    int n = size(A, 1);

    for (int k = 1; k <= s; k++)
    {
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= i - 1; j++)
                B(k, i) = B(k, i) - A(j, i)*B(k, j);
            B(k, i) = B(k, i)/A(i, i);
        }
    }
}

void maxpy(Matrix<double> &A, Matrix<double> &x, Matrix<double> &y)
{
    int m = size(y, 1);
    int l = size(y, 2);
    int n = size(A, 2);

    for (int k = 1; k <= l; k++)
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
                y(i, k) = y(i, k) - A(i, j)*x(j, k);
}