#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "test_cases.h"
#include "containers.h"

template <class T> Vector<T> RandomVector(int length, double nmax)
{
	Vector<T> a(length);
	for (int i = 0; i < length; i++) a[i] = nmax * rand();
	return a;
}

template <class T> Matrix<T> RandomMatrix(int M_rows, int N_cols, int gc, double nmax)
{
	int N_out = N_cols + 2 * gc;
	Matrix<T> m(M_rows, Vector<double>(N_out));
	int N_int = m[0].size() - 2 * gc;

	for (int i = 0; i < M_rows; i++)
		for (int j = gc; j < N_int + gc; j++)
			m[i][j] = nmax * rand();

	return m;
}


#endif