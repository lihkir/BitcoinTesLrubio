#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "test_cases.h"

template <class T> T max(T a, T b)
{
	return a > b ? a : b;
}

template <class T> void PrintingContainer(const std::vector<T>& v, int d)
{
	std::cout << std::fixed;
	std::cout << std::setprecision(d);
	for (auto i = 0; i < v.size(); i++) std::cout << v[i] << "\t";
	std::cout << "\n\n";
}

template <class T> void PrintingContainer(const std::vector<std::vector<T>>& m, int d)
{
	std::cout << std::fixed;
	std::cout << std::setprecision(d);
	for (auto i = 0; i < m.size(); i++)
	{
		for (auto j = 0; j < m[i].size(); j++)
			std::cout << m[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "\n";
}

template <class T> T VectorSum(const std::vector<T>& v)
{
	T s{};
	for (auto i = 0; i != v.size(); i++) s = s + v[i];
	return s;
}

template <class T> std::vector<T> RandomVector(int length, double nmax)
{
	std::vector<T> a(length);
	for (int i = 0; i < length; i++) a[i] = nmax * rand();
	return a;
}

template <class T> std::vector<std::vector<T>> RandomMatrix(int M_rows, int N_cols, double nmax)
{
	struct test_cases* pt_test = get_tests();

	auto gc2 = 2 * pt_test->gc;	auto N_out = N_cols + gc2;
	std::vector<std::vector<T>> m(M_rows, std::vector<double>(N_out));
	auto N_int = m[0].size() - gc2;

	for (auto i = 0; i < M_rows; i++)
		for (auto j = pt_test->gc; j < N_int + pt_test->gc; j++)
			m[i][j] = nmax * rand();

	return m;
}

template <class T> std::vector<T> Col(const std::vector<std::vector<T>>& m, int j)
{
	std::vector<T> c(m.size());
	for (auto i = 0; i < c.size(); i++)
		c[i] = m[i][j];
	return c;
}

template <class T> T DotProduct(const std::vector<T>& u, const std::vector<T>& v)
{
	T s{};
	for (auto i = 0; i < u.size(); i++)
		s += u[i] * v[i];
	return s;
}

template <typename T, typename U> auto max(T x, U y) -> decltype(x > y ? x : y)
{
	return x > y ? x : y;
}

template <class T> std::vector<std::vector<double>> SubMatrix(const std::vector<std::vector<T>>& m)
{
	struct test_cases* pt_test = get_tests();

	auto N_int = (int)m[0].size() - 2 * pt_test->gc;
	std::vector<std::vector<double>> R(m.size(), std::vector<double>(N_int));

	for (auto i = 0; i < R.size(); i++)
	{
		for (auto j = 0; j < R[i].size(); j++)
		{
			R[i][j] = m[i][(size_t)pt_test->gc + j];
		}
	}
	return R;
}

#endif