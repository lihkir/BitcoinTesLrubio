#include "update_utilities.h"

void update_col(const Matrix<double>& u, const Matrix<double>& c, int j)
{
	int N_rows = size(u, 1);
	for (int i = 1; i <= N_rows; i++) u(i, j) = c(i);
}

void copy_block(Block<double>& u, Block<double>& v)
{
	int N_rows = size(u, 1);
	int N_cols = size(u, 2);
	int N_levs = size(u, 3);
	
	for (int i = 1; i <= N_rows; i++)
		for (int j = 1; j <= N_cols; j++)
			for (int k = 1; k <= N_levs; k++)
				u(i, j, k) = v(i, j, k);
}

void update_inside(Matrix<double>& u, Matrix<double>& v, int gc)
{
	int N_rows = size(v, 1); int N_colsv = size(v, 2); int N_colsu = size(u, 2);
	int N_cols, gcl, gcr;

	if (N_colsv < N_colsu)
	{
		N_cols = N_colsv;
		gcl = gc;
		gcr = 0;
	}
	else
	{
		N_cols = N_colsu;
		gcr = gc;
		gcl = 0;
	}
	loop_over(u, v, N_cols, gcl, gcr);
}

void loop_over(Matrix<double> &u, Matrix<double> &v, int N_col, int gcl, int gcr)
{
	int N_rows = size(u, 1);
	for (int i = 1; i <= N_rows; i++)
		for (int j = 1; j <= N_col; j++)
			u(i, j + gcl) = v(i, j + gcr);
}

void get_col(Matrix<double>& c, Matrix<double>& u, int j)
{
	int length = size(c, 1); 
	for (int i = 1; i <= length; i++) c(i) = u(i, j);
}