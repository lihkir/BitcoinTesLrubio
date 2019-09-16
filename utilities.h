#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <iomanip>
#include "test_cases.h"

template <class T> T max(T a, T b)
{
	return a > b ? a : b ;
}

template <class T> void PrintingContainer(const std::vector<T> &v)
{
    std::cout << std::fixed;
	std::cout << std::setprecision(16);
	for (unsigned int i = 0; i < v.size(); i++) std::cout << v[i] << "\t";
	std::cout << "\n\n";
}

template <class T> void PrintingContainer(const std::vector<std::vector<T>> &m)
{
	std::cout << std::fixed;
	std::cout << std::setprecision(16);
	for (unsigned int i = 0; i < m.size(); i++)
	{
		for (unsigned int j = 0; j < m[i].size(); j++)
			std::cout << m[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "\n";
}

template <class T> T VectorSum(const std::vector<T> &v)
{
    T s{};
	for (unsigned int i = 0; i != v.size(); i++) s = s + v[i];
	return s;
}

template <class T> std::vector<T> RandomVector(int length, int nmax)
{
	std::vector<T> a(length);
	for (int i = 0; i < length; i++)
		a[i] = T(rand() % nmax) + 1;
	return a;
}

template <class T> std::vector<std::vector<T>> RandomMatrix(int M_rows, int N_cols, int nmax)
{
	struct test_cases *pt_test = get_tests();

	std::vector<std::vector<T>> m(M_rows, std::vector<double>(N_cols + 2 * pt_test->gc));
	int N_int = m[0].size() - 2 * pt_test->gc;

	for (int i = 0; i < M_rows; i++)
		for (int j = pt_test->gc; j < N_int + pt_test->gc; j++)
			m[i][j] = T(rand() % nmax) + 1;

	return m;
}

template <class T> std::vector<T> Row(const std::vector<std::vector<T>> &m, int j)
{
	std::vector<T> c(m.size());
	for (unsigned int i = 0; i < c.size(); i++)
		c[i] = m[i][j];
	return c;
}

template <class T> T DotProduct(const std::vector<T> &u, const std::vector<T> &v)
{
	T s{};
	for (unsigned int i = 0; i < u.size(); i++)
		s += u[i]*v[i];
	return s;
}

template <typename T, typename U> auto max(T x, U y) -> decltype(x > y ? x : y)
{
  return x > y ? x : y;
}

#endif