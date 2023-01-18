#include "Vector2.h"
#include "Vector3.h"

Vector2::Vector2()
{
	x = 0;
	y = 0;
}

Vector2::Vector2(double _x, double _y)
{
	x = _x;
	y = _y;
}

Vector2::~Vector2()
{
}

Vector3 Vector2::Homogeneous() const
{
	return Vector3(x, y, 1.0);
}

void Vector2::Print() const
{
	std::cout << "Vector2(" << x << "," << y << ")" << std::endl;
}
