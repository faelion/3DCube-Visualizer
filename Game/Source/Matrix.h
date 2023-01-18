#pragma once

#include "Array.h"

template<class T, int R, int C>
class Matrix
{
private:
	Array<T, C> m_Matrix[R];

	double det_impl(const T* data, int rows, int columns) const;
public:
	Matrix();
	Matrix(T* arr);

	int& at(int r, int c){ return m_Matrix[r][c]; }
	Array<T, C>& operator[](int index) { return m_Matrix[index]; }
	const Array<T, C>& operator[](int index) const { return m_Matrix[index]; }

	double Trace() const;
	T Determinant() const;
	Matrix<T, R, C> Transpose() const;
	Matrix<T, R, C> Inverse() const;
	Matrix<T, R, C> Adjoint() const;

	Matrix<T, R, C> Identity() const;

	int RowSize() { return R; }
	int ColumnSize() { return C; }
	int Size() { return R * C; }
	void FillMatrix(T* matrix = 0);

	void Print() const;
};

template<class T, int R, int C>
inline Matrix<T, R, C> operator*(const Matrix<T, R, C>& m1, const Matrix<T, R, C>& m2) {
	Matrix<T, R, C> mult;

	for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < C; j++)
		{
			for (int k = 0; k < R; k++)
			{
				mult[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}

	return mult;
}

template<class T, int R, int C>
inline double Matrix<T, R, C>::det_impl(const T* data, int rows, int columns) const
{
	// If size of rows = 2, calculate determinant of the matrix
	if (rows == 2) {
		return data[0] * data[3] - data[1] * data[2];
	}

	double determinant = 0;
	int sign = 1;

	int subW = rows - 1;
	int subH = columns - 1;

	T* subMatrix = new T[subW * subH];

	int i = 0;
	for (int j = 0; j < columns; j++) {
		int p = 0;
		// Create a submatrix
		for (int ii = 0; ii < rows; ii++) {
			int q = 0;
			if (ii == i)
				continue;

			for (int jj = 0; jj < columns; jj++) {
				if (jj == j)
					continue;

				subMatrix[p * subW + q] = data[ii * rows + jj];
				q++;
			}

			p++;
		}

		// Sum determinants
		determinant += sign * data[i * rows + j] * det_impl(subMatrix, subW, subH);
		sign *= -1;
	}

	delete[] subMatrix;

	return determinant;
}

template<class T, int R, int C>
inline Matrix<T, R, C>::Matrix()
{
	FillMatrix(0);
}

template<class T, int R, int C>
inline Matrix<T, R, C>::Matrix(T * arr)
{
	FillMatrix(arr);
}

template<class T, int R, int C>
inline double Matrix<T, R, C>::Trace() const
{
	double trace = 0.0;

	for (int i = 0; i < R; i++) {
		trace += m_Matrix[i][i];
	}

	return trace;
}

template<class T, int R, int C>
inline T Matrix<T, R, C>::Determinant() const
{
	return det_impl(&m_Matrix[0][0], R, C);
}

template<class T, int R, int C>
inline Matrix<T, R, C> Matrix<T, R, C>::Transpose() const
{
	Matrix<T, R, C> transpose;

	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			transpose[i][j] = m_Matrix[j][i];
		}
	}

	return transpose;
}

template<class T, int R, int C>
inline Matrix<T, R, C> Matrix<T, R, C>::Inverse() const
{
	Matrix<T, R, C> inv;

	Matrix<T, R, C> adjoint = Adjoint();
	double det = Determinant();

	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			inv[i][j] = adjoint[i][j] / det;
		}
	}

	return inv;
}

template<class T, int R, int C>
inline Matrix<T, R, C> Matrix<T, R, C>::Adjoint() const
{
	Matrix<T, R, C> adjoint;

	Matrix<T, R - 1, C - 1> subAdj;

	int sign = 1;

	// Find cofactor of each value in matrix
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			int p=0;
			// Create a submatrix
			for (int ii = 0; ii < R; ii++) {
				int q = 0;
				if (ii == i)
					continue;

				for (int jj = 0; jj < C; jj++) {
					if (jj == j)
						continue;

					subAdj[p][q] = m_Matrix[ii][jj];
					q++;
				}

				p++;
			}

			// Calculate cofactor
			adjoint[i][j] = sign * subAdj.Determinant();
			sign *= -1;
		}
	}

	return adjoint.Transpose();
}

template<class T, int R, int C>
inline Matrix<T, R, C> Matrix<T, R, C>::Identity() const
{
	Matrix<T, R, C> identity;

	for (int i = 0; i < R; i++) {
		identity[i][i] = 1;
	}

	return identity;
}

template<class T, int R, int C>
inline void Matrix<T, R, C>::FillMatrix(T * matrix)
{
	if (matrix) {
		for (int i = 0; i < R; i++) {
			m_Matrix[i].FillArray(&matrix[i * R]);
		}
	}
	else {
		for (int i = 0; i < R; i++) {
			m_Matrix[i].FillArray(0);
		}
	}
}

template<class T, int R, int C>
inline void Matrix<T, R, C>::Print() const
{
	for (int i = 0; i < R; i++) {
		m_Matrix[i].Print();
	}
}
