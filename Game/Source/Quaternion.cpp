#include "Quaternion.h"
#include "Vector3.h"
#include <math.h>


Quaternion::Quaternion()
{
}

Quaternion::Quaternion(double _x, double _y, double _z, double _w)
{
	x = _x;
	y = _y;
	z = _z;
	w = _w;
}


Quaternion::~Quaternion()
{
}

void Quaternion::Normalize()
{
	double n = Magnitude();

	x /= n;
	y /= n;
	z /= n;
	w /= n;
}

Quaternion Quaternion::Normalized() const
{
	Quaternion q;

	double n = Magnitude();

	q.x /= n;
	q.y /= n;
	q.z /= n;
	q.w /= n;

	return q;
}

void Quaternion::Scale(double s)
{
	x *= s;
	y *= s;
	z *= s;
	w *= s;
}

Quaternion Quaternion::Conjugate() const
{
	Quaternion quat;

	quat.x = -x;
	quat.y = -y;
	quat.z = -z;
	quat.w = w;

	return quat;
}

Quaternion Quaternion::Inverse() const
{
	Quaternion q(*this);

	double m = Magnitude();

	if (m != 0) {

		Quaternion conjugate = Conjugate();

		q.x = conjugate.x / (m * m);
		q.y = conjugate.y / (m * m);
		q.z = conjugate.z / (m * m);
		q.w = conjugate.w / (m * m);
	} // else no inverse, return same quaternion

	return q;
}

void Quaternion::SetRotation(double _x, double _y, double _z, double angle)
{
	Vector3 v(_x, _y, _z);
	v.Normalize();

	double c = cos(angle / 2);
	double s = sin(angle / 2);

	x = v.x * s;
	y = v.y * s;
	z = v.z * s;
	w = c;
}

double Quaternion::Magnitude() const
{
	return sqrt(x * x + y * y + z * z + w * w);
}

void Quaternion::Print() const
{
	std::cout << "Quaternion(" << x << "," << y << "," << z << "," << w << ")" << std::endl;
}

Quaternion operator*(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q;

	q.x = q1.w * q2.x + q1.y * q2.z - q1.z * q2.y + q1.x * q2.w;
	q.y = q1.w * q2.y + q1.y * q2.w + q1.z * q2.x - q1.x * q2.z;
	q.z = q1.w * q2.z - q1.y * q2.x + q1.z * q2.w + q1.x * q2.y;
	q.w = q1.w * q2.w - q1.y * q2.y - q1.z * q2.z - q1.x * q2.x;

	return q;
}

Quaternion operator+(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q;

	q.x = q1.x + q2.x;
	q.y = q1.y + q2.y;
	q.z = q1.z + q2.z;
	q.w = q1.w + q2.w;

	return q;
}
