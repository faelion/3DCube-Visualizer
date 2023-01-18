#pragma once

class Vector3;

#include <iostream>

class Quaternion
{
public:
	Quaternion();
	Quaternion(double _x, double _y, double _z, double _w);
	~Quaternion();

	void Normalize();
	Quaternion Normalized() const;

	void Scale(double s);

	Quaternion Conjugate() const;
	Quaternion Inverse() const;

	void SetRotation(double _x, double _y, double _z, double angle);

	double Magnitude() const;

	void Print() const;

	double x, y, z, w;
};

Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
Quaternion operator+(const Quaternion& q1, const Quaternion& q2);