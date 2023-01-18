#pragma once

#include <iostream>

class Vector3;

class Vector2
{
public:
	Vector2();
	Vector2(double _x, double _y);
	~Vector2();

	double x;
	double y;

	Vector3 Homogeneous() const;

	void Print() const;
};

