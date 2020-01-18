#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <string.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

/*****************************************/
/** Vector class "Vector" **/

#include <cassert>

template <class T> class Matrix;
template <class T> class Block;

template <class T> class Vector
{
	private:
		int size;
		T* data;
		
	public:
		Vector();
		Vector(int);
		Vector(int, T);
		Vector(const Vector<T>&);
		~Vector();

		T &operator [] (int) const;
		T &operator () (int) const;
		T operator | (const Vector<T>&);
		Vector<T> operator *= (T);
		Vector<T> operator * (T) const;
		Vector<T> operator + (const Vector<T>&);

		int length() const;
		void resize(int);
		void resize(int, T);

		friend class Matrix<T>;
		
		template <class U> friend void print1D(const Vector<U>& v, int dp, std::string name);
		template <class U> friend U sum1D(const Vector<U>& v);
		template <class U> friend U dot_product(const Vector<U>& u, const Vector<U>& v);
		template <class U> friend U square_sum(const Vector<U>& u);
};

template <class T> Vector<T>::Vector() : size(0), data(NULL) {}

template <class T> Vector<T>::Vector(int n) : size(n), data(new T[n]) { assert(data != NULL); }

template <class T> Vector<T>::Vector(int n, T value) : size(n), data(new T[n])
{
	assert(data != NULL);
	for (int i = 0; i < n; i++) data[i] = value;
}

template <class T> Vector<T>::Vector(const Vector<T>& v) : size(v.size), data(new T[v.size])
{
	assert(data != NULL);
	for (int i = 0; i < v.size; i++) data[i] = v.data[i];
}

template <class T> Vector<T>::~Vector() { delete[] data; }

template <class T> T &Vector<T>::operator[](int i) const
{
	assert(i >= 0 && i < size);
	return data[i];
}

template <class T> T &Vector<T>::operator()(int i) const
{
	assert(i >= 1 && i =< size);
	return data[i - 1];
}

template <class T> T Vector<T>::operator | (const Vector<T>& v)
{
	assert(size == v.size);
	T result{};

	for (int i = 0; i < size; i++) result = result + data[i] * v.data[i];
	return result;
}

template <class T> Vector<T> Vector<T>::operator *= (T c)
{
	for (int i = 0; i < size; i++) data[i] *= c;
	return *this;
}

template <class T> Vector<T> Vector<T>:: operator * (T c) const
{
	Vector<T> result(*this);
	return result *= c;
}

template <class T> int Vector<T>::length() const { return size; }

template <class T> void print1D(const Vector<T> &v, int dp, std::string name)
{
	std::cout << name << " = \n\n";
	std::cout << std::fixed;
	std::cout << std::setprecision(dp);
	for (int i = 0; i < v.length(); i++) std::cout << v[i] << "\n";
	std::cout << "\n\n";
}

template <class T> T sum1D(const Vector<T>& v)
{
	T s{};
	for (auto i = 0; i != v.length(); i++) s = s + v[i];
	return s;
}

template <class T> void Vector<T>::resize(int length)
{
	int i;
	T zero(0);
	T* newData = new T[length]; assert(newData != NULL);

	if (length <= size)
		for (i = 0; i < length; i++) newData[i] = data[i];
	else
	{
		for (i = 0; i < size; i++) newData[i] = data[i];
		for (i = size; i < length; i++) newData[i] = zero;
	}

	delete[] data;
	size = length;
	data = newData;
}

template <class T> void Vector<T>::resize(int length, T value)
{
	int i;
	T* newData = new T[length]; assert(newData != NULL);

	if (length <= size)
		for (i = 0; i < length; i++) newData[i] = data[i];
	else
	{
		for (i = 0; i < size; i++) newData[i] = data[i];
		for (i = size; i < length; i++) newData[i] = value;
	}
	delete[] data;
	size = length;
	data = newData;
}

template <class T> Vector<T> Vector<T>::operator + (const Vector<T>& v)
{
	assert(size == v.size);
	Vector<T>* result = new Vector<T>(size);
	for (int i = 0; i < size; i++) result->data[i] = v.data[i] + data[i];
	return *result;
}

template <class T> T dot_product(const Vector<T>& u, const Vector<T>& v)
{
	T s{}; int N = u.length();
	for (int i = 0; i < N; i++) s += u[i] * v[i];
	return s;
}

template <class T> T square_sum(const Vector<T>& u)
{
	T s{}; int N = u.length();
	for (int i = 0; i < N; i++) s += u[i] * u[i];
	return s;
}

template <class T> Vector<T> sub_vector(const Matrix<T>& m, int col)
{
	Vector<T> R(m.rows());
	for (unsigned int i = 0; i < R.size(); i++) R[i] = m[i][col];
	return R;
}

/*****************************************/
/** Matrix class "Matrix" **/

template <class T> class Matrix
{
	private:
		int rowNum, colNum;
		Vector<T>* mat;
	public:
		Matrix();
		Matrix(int, int);
		Matrix(int, int, T);
		Matrix(const Vector<T>&);
		Matrix(const Matrix<T>&);
		~Matrix();

		void resize(int, int);
		void resize(int, int, T);
		
		const Matrix<T>& operator = (const Matrix<T>&);
		Vector<T> &operator [] (int) const;
		T &operator () (int, int) const;
		T &operator () (int) const;
		Matrix<T> operator * (const Matrix<T>&) const;
		Matrix<T> operator *= (T) const;
		Matrix<T> operator * (T) const;
		Matrix<T> operator + (const Matrix<T>&) const;

		Matrix<T> transpose() const;
		inline int rows() const;
		inline int cols() const;

		friend class Block<T>;

		template <class U> friend void print2D(const Matrix<U> &v, int dp, std::string name);
		template <class U> friend Matrix<U> sub_matrix(const Matrix<U>& m, int gc);
		template <class U> friend U sum1D(const Matrix<U>& v);
		template <class U> friend U dot_product(const Matrix<U>& u, const Matrix<U>& v);
		template <class U> friend inline int size(const Matrix<U>& u, int idx);
};

template <class T> Matrix<T>::Matrix() : rowNum(0), colNum(0), mat(NULL) {}

template <class T> Matrix<T>::Matrix(int r, int c) : rowNum(r), colNum(c), mat(new Vector<T>[r])
{
	assert(mat != NULL);
	for (int i = 0; i < r; i++) mat[i].resize(c);
}

template <class T> Matrix<T>::Matrix(int r, int c, T value) : rowNum(r), colNum(c), mat(new Vector<T>[r])
{
	assert(mat != NULL);
	for (int i = 0; i < r; i++) mat[i].resize(c, value);
}

template <class T> Matrix<T>::Matrix(const Vector<T> &v) : rowNum(v.length()), colNum(1), mat(new Vector<T>[rowNum])
{
	assert(mat != NULL);
	for (int i = 0; i < rowNum; i++) mat[i].resize(1, v[i]);
}

template <class T> Matrix<T>::Matrix(const Matrix<T> &m) : rowNum(m.rowNum), colNum(m.colNum), mat(new Vector<T>[m.rowNum])
{
	assert(mat != NULL);
	for (int i = 0; i < m.rowNum; i++) mat[i] = m.mat[i];
}

template <class T> Matrix<T>::~Matrix() { delete[] mat; }

template <class T> const Matrix<T>& Matrix<T>::operator = (const Matrix<T>& m)
{
	if (this == &m) return *this;
	delete[] mat;
	rowNum = m.rowNum; colNum = m.colNum;
	mat = new Vector<T>[m.rowNum]; assert(mat != NULL);
	for (int i = 0; i < m.rowNum; i++) mat[i] = m.mat[i];
	return *this;
}

template <class T> Vector<T> &Matrix<T>::operator [] (int index) const
{
	assert(index >= 0 && index < rowNum);
	return mat[index];
}

template <class T> T &Matrix<T>::operator () (int index) const
{
	if (colNum == 1)
		return mat[index - 1][0];
	else
		return mat[0][index -1];
}

template <class T> T &Matrix<T>::operator () (int i, int j) const
{
	assert((i >= 1 && i <= rowNum) && (j >= 1 && j <= colNum));
	return mat[i - 1][j - 1];
}

template <class T> Matrix<T> Matrix<T>::operator * (const Matrix<T>& m) const
{
	assert(colNum == m.rowNum);
	Matrix<T> *result = new Matrix<T>(rowNum, m.colNum);

	for (int i = 0; i < rowNum; i++)
		for (int j = 0; j < m.colNum; j++)
			result->mat[i][j] = mat[i] | m(j);
	
	return *result;
}

template <class T> Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T>* result = new Matrix<T>(colNum, rowNum);

	for (int i = 0; i < rowNum; i++)
		for (int j = 0; j < colNum; j++)
			result->mat[j][i] = mat[i][j];

	return *result;
}

template <class T> int Matrix<T>::rows() const { return rowNum; }

template <class T> int Matrix<T>::cols() const { return colNum; }

template <class T> void print2D(const Matrix<T> &m, int dp, std::string name)
{
	std::cout << name << " = \n\n";
	std::cout << std::fixed;
	std::cout << std::setprecision(dp);
	for (int i = 0; i < m.rows(); i++)
	{
		for (int j = 0; j < m.cols(); j++)
			std::cout << m[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "\n";
}

template <class T> Matrix<T> Matrix<T>::operator *= (T c) const
{
	for (int i = 0; i < rowNum; i++) mat[i] *= c;
	return *this;
}

template <class T> Matrix<T> Matrix<T>::operator * (T value) const
{
	Matrix<T>* result = new Matrix<T>(rowNum, colNum);
	
	for (int i = 0; i < rowNum; i++)
		for (int j = 0; j < colNum; j++)
			result->mat[i][j] = mat[i][j]*value;
	
	return *result;
}

template <class T> Matrix<T> Matrix<T>::operator + (const Matrix<T>& m) const
{
	assert(rowNum == m.rowNum && colNum == m.colNum);
	Matrix<T>* result = new Matrix<T>(rowNum, colNum);

	for (int i = 0; i < rowNum; i++)
		for (int j = 0; j < colNum; j++)
			result->mat[i][j] = mat[i][j] + m.mat[i][j];

	return *result;
}

template<class T> void Matrix<T>::resize(int r, int c)
{
	int i;
	Vector<T> *newMatrix = new Vector<T>[r]; assert(newMatrix != NULL);
	
	if (r <= rowNum)
	{
		for (int i = 0; i < r; i++)
		{
			(mat + i)->resize(c);
			newMatrix[i] = mat[i];
		}
	}
	else
	{
		for (int i = 0; i < rowNum; i++)
		{
			(mat + i)->resize(c);
			newMatrix[i] = mat[i];
		}
		for (int i = rowNum; i < r; i++) newMatrix[i].resize(c);		
	}
	delete [] mat;
	rowNum = r; colNum = c;
	mat = newMatrix;
}

template<class T> void Matrix<T>::resize(int r, int c, T value)
{
	int i;
	Vector<T> *newMatrix = new Vector<T>[r]; assert(newMatrix != NULL);
	
	if (r <= rowNum)
	{
		for (int i = 0; i < r; i++)
		{
			(mat + i)->resize(c, value);
			newMatrix[i] = mat[i];
		}
	}
	else
	{
		for (int i = 0; i < rowNum; i++)
		{
			(mat + i)->resize(c, value);
			newMatrix[i] = mat[i];
		}
		for (int i = rowNum; i < r; i++) newMatrix[i].resize(c, value);		
	}
	delete [] mat;
	rowNum = r; colNum = c;
	mat = newMatrix;
}

template <class T> Matrix<T> sub_matrix(const Matrix<T>& m, int gc)
{
	int N_int = m.cols() - 2 * gc, M_rows = m.rows(); 
	Matrix<T> R(M_rows, N_int);

	for (int i = 0; i < M_rows; i++)
		for (int j = 0; j < N_int; j++)
			R[i][j] = m[i][gc + j];

	return R;
}

template <class T> T sum1D(const Matrix<T>& v)
{
	T s{};
	if (v.cols() == 1)
		for (auto i = 0; i != v.rows(); i++) s = s + v[i][0];
	else
		for (auto i = 0; i != v.cols(); i++) s = s + v[0][i];
	return s;
}

template <class T> T dot_product(const Matrix<T>& u, const Matrix<T>& v)
{
	T s{}; 
	if (u.cols() == 1)
		for (int i = 0; i < u.rows(); i++) s += u[i][0] * v[i][0];
	else
		for (int i = 0; i < u.cols(); i++) s += u[0][i] * v[0][i];
	return s;
}

template <class U> int size(const Matrix<U>& u, int idx) { return idx == 1 ? u.rows() : u.cols(); }

/*****************************************/
/** Block class "Block" **/

template <class T> class Block
{
	private:
		int rowNum, colNum, levNum;
		Matrix<T>* block;
	public:
		Block(int = 0, int = 0, int = 0);
		Block(int, int, int, T);
		Block(const Block<T>&);
		~Block();

		const Block<T>& operator = (const Block<T>&);
		Matrix<T> &operator [] (int) const;
		T &operator () (int, int, int) const;

		inline int rows() const;
		inline int cols() const;
		inline int levs() const;

		template <class U> friend inline int size(const Block<U>& u, int idx);
};

template <class T> Block<T>::Block(int r, int c, int v) : rowNum(r), colNum(c), levNum(v), block(new Matrix<T>[r])
{
	assert(block != NULL);
	for (int i = 0; i < r; i++) (block + i)->resize(c, v);
}

template <class T> Block<T>::Block(int r, int c, int v, T num) : rowNum(r), colNum(c), levNum(v), block(new Matrix<T>[r])
{
	assert(block != NULL);
	for (int i = 0; i < r; i++) (block + i)->resize(c, v, num);
}

template <class T> Block<T>::Block(const Block<T> &b) : rowNum(b.rowNum), colNum(b.colNum), levNum(b.levNum), block(new Matrix<T>[b.rowNum])
{
	assert(block != NULL);
	for (int i = 0; i < b.rowNum; i++) block[i] = b.block[i];
}

template <class T> Block<T>::~Block() { delete[] block; }

template <class T> Matrix<T> &Block<T>::operator [] (int index) const
{
	assert(index >= 0 && index < rowNum);
	return block[index];
}

template <class T> T &Block<T>::operator () (int b, int i, int j) const
{
	assert((b >= 1 && b <= rowNum) && (i >= 1 && i <= colNum) && (j >= 1 && j <= levNum));
	return block[b - 1][i - 1][j - 1];
}

template <class T> int Block<T>::rows() const { return rowNum; }

template <class T> int Block<T>::cols() const { return colNum; }

template <class T> int Block<T>::levs() const { return levNum; }

template <class U> int size(const Block<U>& u, int idx) { return idx == 1 ? u.rows() : idx == 2 ? u.cols() : u.levs(); }

template <class T> const Block<T> &Block<T>::operator = (const Block<T>& m)
{
	if (this == &m) return *this;
	delete[] block;
	rowNum = m.rowNum; colNum = m.colNum; levNum = m.levNum;
	block = new Matrix<T>[m.rowNum]; assert(mat != NULL);
	for (int i = 0; i < m.rowNum; i++) block[i] = m.block[i];
	return *this;
}

/*****************************************/
/** Usefull Functions **/

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

#endif