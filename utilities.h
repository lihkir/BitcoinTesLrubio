#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "test_cases.h"

template <class T> inline T max(T a, T b)
{
	return a > b ? a : b;
}

template <class T> inline T min(T a, T b)
{
	return a < b ? a : b;
}

template <typename T> inline int sgn(T val) 
{
	return (T(0) < val) - (val < T(0));
}

template <class T> void PrintingContainer(const std::vector<T>& v, int d)
{
	std::cout << std::fixed;
	std::cout << std::setprecision(d);
	for (unsigned int i = 0; i < v.size(); i++) std::cout << v[i] << "\t";
	std::cout << "\n\n";
}

template <class T> void PrintingContainer(const std::vector<std::vector<T>>& m, int d)
{
	std::cout << std::fixed;
	std::cout << std::setprecision(d);
	for (unsigned int i = 0; i < m.size(); i++)
	{
		for (unsigned int j = 0; j < m[i].size(); j++)
			std::cout << m[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "\n";
}

template <class T> T VectorSum(const std::vector<T>& v)
{
	T s{};
	for (unsigned int i = 0; i != v.size(); i++) s = s + v[i];
	return s;
}

template <class T> std::vector<T> RandomVector(int length, double nmax)
{
	std::vector<T> a(length);
	for (int i = 0; i < length; i++) a[i] = nmax * rand();
	return a;
}

template <class T> std::vector<std::vector<T>> RandomMatrix(int M_rows, int N_cols, int gc, double nmax)
{
	int N_out = N_cols + 2 * gc;
	std::vector<std::vector<T>> m(M_rows, std::vector<double>(N_out));
	int N_int = m[0].size() - 2 * gc;

	for (int i = 0; i < M_rows; i++)
		for (int j = gc; j < N_int + gc; j++)
			m[i][j] = nmax * rand();

	return m;
}

template <class T> std::vector<T> Col(const std::vector<std::vector<T>>& m, int j)
{
	std::vector<T> c(m.size());
	for (unsigned int i = 0; i < c.size(); i++)
		c[i] = m[i][j];
	return c;
}

template <class T> T DotProduct(const std::vector<T>& u, const std::vector<T>& v)
{
	T s{};
	for (unsigned int i = 0; i < u.size(); i++)
		s += u[i] * v[i];
	return s;
}

template <class T> T SquareSum(const std::vector<T>& u)
{
	T s{};
	for (unsigned int i = 0; i < u.size(); i++)
		s += u[i] * u[i];
	return s;
}

template <typename T, typename U> auto max(T x, U y) -> decltype(x > y ? x : y)
{
	return x > y ? x : y;
}

template <class T> std::vector<std::vector<double>> SubMatrix(const std::vector<std::vector<T>>& m, int gc)
{
	int N_int = (int)m[0].size() - 2 * gc;
	std::vector<std::vector<double>> R(m.size(), std::vector<double>(N_int));

	for (unsigned int i = 0; i < R.size(); i++)
		for (unsigned int j = 0; j < R[i].size(); j++)
			R[i][j] = m[i][gc + j];

	return R;
}

template <class T> std::vector<T> SubVector(const std::vector<std::vector<T>>& m, int col)
{
	std::vector<T> R(m.size());
	for (unsigned int i = 0; i < R.size(); i++) R[i] = m[i][col];
	return R;
}

template <class T> std::vector<std::vector<T>> SubCol(const std::vector<std::vector<T>>& m, int col)
{
	std::vector<std::vector<T>> R(m.size(), std::vector<T>(1));
	for (unsigned int i = 0; i < R.size(); i++) R[i][0] = m[i][col];
	return R;
}

#endif