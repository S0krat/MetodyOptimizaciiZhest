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

	double operator*(const Point p) {
		return x * p.x + y * p.y;
	}
};

class Matrix {
public:
	Point a;
	Point b;

	Matrix() {
		a = Point();
		b = Point();
	} 

	Matrix(Point p1, Point p2) {
		a = p1;
		b = p2;
	}

	double det() {
		return a.x * b.y - a.y * b.x;
	}

	void inverse() {
		double det = this->det();
		if (det == 0) return;
		double temp = a.x;
		a.x = b.y / det;
		b.y = temp / det;
		a.y /= -det;
		b.x /= -det;
	}

	Point operator*(const Point p) {
		return Point(a.x * p.x + a.y * p.y, b.x * p.x + b.y * p.y);
	}
};