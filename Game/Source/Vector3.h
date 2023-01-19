#pragma once

class Quaternion;
class Vector4;

#include <iostream>

class Vector3
{
public:
	Vector3();
	Vector3(double _x, double _y, double _z);
	~Vector3();

	double x, y, z;

	void Normalize();
	void Rotate(const Quaternion& q);
	void Rotate(const Vector3& axis, double angle);
	void Rotate(float x, float y, float z);
	Vector4 Homogeneous() const;

	bool IsNormalized();

	double Magnitude();
	void Print() const;
};

double Dot(const Vector3& v1, const Vector3& v2);
Vector3 Cross(const Vector3& v1, const Vector3& v2);

Vector3 operator+(const Vector3& v1, const Vector3& v2);

Vector3 operator+(const Vector3& v, double s);
Vector3 operator*(double s, const Vector3& v);
Vector3 operator*(const Vector3& v, double s);
Vector3 operator/(const Vector3& v, double s);