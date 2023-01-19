#include "Vector3.h"
#include "Quaternion.h"
#include "Vector4.h"
#include "Math.h"

#include <iostream>

Vector3::Vector3()
{
}

Vector3::Vector3(double _x, double _y, double _z)
{
	x = _x;
	y = _y;
	z = _z;
}


Vector3::~Vector3()
{
}

void Vector3::Normalize()
{
	double m = Magnitude();

	if (m != 0) {
		x /= m;
		y /= m;
		z /= m;
	}
}

void Vector3::Rotate(const Quaternion & q)
{
	// By vector calculations
	/*Vector3 u(q.x, q.y, q.z);
	Vector3 v(x, y, z);

	float s = q.w;

	Vector3 result = 2.0f * Dot(u, v) * u
		+ (s*s - Dot(u, u)) * v
		+ 2.0f * s * Cross(u, v);

	x = result.x;
	y = result.y;
	z = result.z;*/

	// By quaternion multiplications
	Quaternion qv(x, y, z, 0);
	Quaternion qinv = q.Inverse();

	Quaternion qvprime = q * qv * qinv;
	// w = 0
	x = qvprime.x;
	y = qvprime.y;
	z = qvprime.z;
}

void Vector3::Rotate(const Vector3& axis, double angle)
{
	x = x * cos(angle) + sin(angle) * (x * axis.x) + (1 - cos(angle)) * (x * axis.x) * axis.x;
	y = y * cos(angle) + sin(angle) * (y * axis.y) + (1 - cos(angle)) * (y * axis.y) * axis.y;
	z = z * cos(angle) + sin(angle) * (z * axis.z) + (1 - cos(angle)) * (z * axis.z) * axis.z;
}

void Vector3::Rotate(float a, float b, float c)
{
	float rad = 0;
	rad = a;
	y = cos(rad) * y - sin(rad) * z;
	z = sin(rad) * y + cos(rad) * z;

	rad = b;
	x = cos(rad) * x + sin(rad) * z;
	z = -sin(rad) * x + cos(rad) * z;

	rad = c;
	x = cos(rad) * x - sin(rad) * y;
	y = sin(rad) * x + cos(rad) * y;

}

Vector4 Vector3::Homogeneous() const
{
	return Vector4(x, y, z, 1.0);
}

bool Vector3::IsNormalized()
{
	return Magnitude() == 1 ? true : false;
}

double Vector3::Magnitude()
{
	return sqrt(x * x + y * y + z * z);
}

void Vector3::Print() const
{
	std::cout << "Vector3(" << x << "," << y << "," << z << ")" << std::endl;
}

double Dot(const Vector3 & v1, const Vector3 & v2)
{
	double product = 0.0;

	product += v1.x * v2.x;
	product += v1.y * v2.y;
	product += v1.z * v2.z;

	return product;
}

Vector3 Cross(const Vector3 & v1, const Vector3 & v2)
{
	Vector3 v;

	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = -(v1.x * v2.z - v1.z * v2.x);
	v.z = v1.x * v2.y - v1.y * v2.x;

	return v;
}

Vector3 operator+(const Vector3 & v1, const Vector3 & v2)
{
	Vector3 v;

	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	v.z = v1.z + v2.z;

	return v;
}

Vector3 operator+(const Vector3 & v, double s)
{
	Vector3 vec;

	vec.x = v.x + s;
	vec.y = v.y + s;
	vec.z = v.z + s;

	return vec;
}

Vector3 operator*(double s, const Vector3 & v)
{
	Vector3 vec;

	vec.x = v.x * s;
	vec.y = v.y * s;
	vec.z = v.z * s;

	return vec;
}

Vector3 operator*(const Vector3 & v, double s)
{
	Vector3 vec;

	vec.x = v.x * s;
	vec.y = v.y * s;
	vec.z = v.z * s;

	return vec;
}

Vector3 operator/(const Vector3 & v, double s)
{
	Vector3 vec;

	vec.x = v.x / s;
	vec.y = v.y / s;
	vec.z = v.z / s;

	return vec;
}
