#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include "globals.h"

template <class T> T max(T a, T b)
{
	return a > b ? a : b ;
}

template <class T> void PrintingContainer(const std::vector<T> &v)
{
    std::cout << std::fixed;
	for (unsigned int i = 0; i < v.size(); i++) std::cout << v[i] << "\t";
	std::cout << "\n\n";
}

template <class T> void PrintingContainer(const std::vector<std::vector<T>> &m)
{
	std::cout << std::fixed;
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
	global_data *pt_data;
	std::vector<std::vector<T>> m(M_rows, std::vector<double>(N_cols + 2 * pt_data->gc));
	int N_int = m[0].size() - 2 * pt_data->gc;

	for (int i = 0; i < M_rows; i++)
		for (int j = pt_data->gc; j < N_int + pt_data->gc; j++)
			m[i][j] = T(rand() % nmax) + 1;

	return m;
}

#endif