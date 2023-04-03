#pragma once
#include <math.h>

class Point {
public:
	double x;
	double y;

	Point() {
		x = 0;
		y = 0;
	}

	Point(double x, double y) {
		this->x = x;
		this->y = y;
	}

	double distTo(Point p) {
		return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
	}

	double norm() {
		return sqrt(x * x + y * y);
	}

	Point operator+(const Point p) {
		return Point(x + p.x, y + p.y);
	}

	Point operator-(const Point p) {
		return Point(x - p.x, y - p.y);
	}

	Point operator*(const double a) {
		return Point(a * x, a * y);
	}

	Point operator/(const double a) {
		return Point(x / a, y / a);
	}
};

