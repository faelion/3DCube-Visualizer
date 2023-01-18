#pragma once

#include <iostream>

template<class T, int size>
class Array {
private:
	T m_Array[size];
public:
	Array();
	Array(T* arr);

	T& operator[](int index) { return m_Array[index]; }
	const T& operator[](int index) const { return m_Array[index]; }

	Array<T, size> operator*(T iterm);
	Array<T, size> operator/(T iterm);
	Array<T, size> operator-(const Array<T, size>& arr);

	void operator-=(const Array<T, size>& arr);

	int Size() { return size; }

	void Multiply(T iterm);
	void Divide(T iterm);
	void FillArray(T* arr = 0);

	void Print() const;
};

template<class T, int size>
inline Array<T, size>::Array()
{
	FillArray(0);
}

template<class T, int size>
inline Array<T, size>::Array(T * arr)
{
	FillArray(arr);
}

template<class T, int size>
inline Array<T, size> Array<T, size>::operator*(T iterm)
{
	Array<T, size> arr(*this);
	arr.Multiply(iterm);
	return arr;
}

template<class T, int size>
inline Array<T, size> Array<T, size>::operator/(T iterm)
{
	Array<T, size> arr(*this);
	arr.Divide(iterm);
	return arr;
}

template<class T, int size>
inline Array<T, size> Array<T, size>::operator-(const Array<T, size>& arr)
{
	Array<T, size> res(*this);
	for (int i = 0; i < size; i++) {
		res[i] -= arr[i];
	}
	return res;
}

template<class T, int size>
inline void Array<T, size>::operator-=(const Array<T, size>& arr)
{
	for (int i = 0; i < size; i++) {
		m_Array[i] -= arr[i];
	}
}

template<class T, int size>
inline void Array<T, size>::Multiply(T iterm)
{
	for (int i = 0; i < size; i++) {
		m_Array[i] *= iterm;
	}
}

template<class T, int size>
inline void Array<T, size>::Divide(T iterm)
{
	for (int i = 0; i < size; i++) {
		m_Array[i] /= iterm;
	}
}

template<class T, int size>
inline void Array<T, size>::FillArray(T* arr)
{
	if (arr) {
		for (int i = 0; i < size; i++) {
			m_Array[i] = arr[i];
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			m_Array[i] = 0;
		}
	}
}

template<class T, int size>
inline void Array<T, size>::Print() const
{
	std::cout << "[ ";

	for (int i = 0; i < size; i++) {
		std::cout << m_Array[i] << " ";
	}

	std::cout << "]" << std::endl;
}
