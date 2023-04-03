#include "Point.h"
#include <iostream>

using namespace std;

double function(Point p) {
	return p.x * p.x + 2 * p.y * p.y + sin(5 * p.x + 3 * p.y) + 3 * p.x + 2 * p.y;
}

Point gradient(Point p) {
	return Point(2 * p.x + 5 * cos(5 * p.x + 3 * p.y) + 3, 4 * p.y + 3 * cos(5 * p.x + 3 * p.y) + 2);
}

Point GoldenRatio(Point left, Point right, double eps) {
	const double phi = (3 - sqrt(5)) / 2;
	Point y, z;
	double Fy, Fz;
	int funcCalls = 0;

	y = right - left;
	y = y * phi;
	y = y + left;

	Fy = function(y);
	funcCalls++;

	z = right + (left - y);
	Fz = function(z);
	funcCalls++;

	while (true) {
		if (Fy < Fz) {
			right = z;
			if (right.distTo(left) < eps) break;
			z = y;
			Fz = Fy;
			y = right + (left - y);
			Fy = function(y);
			funcCalls++;
		}
		else {
			left = y;
			if (right.distTo(left) < eps) break;
			y = z;
			Fy = Fz;
			z = left + (right - z);
			Fz = function(z);
			funcCalls++;
		}
	}

	Point ans = (right + left) / 2.;
	return ans;
}

Point Dichotomy(Point left, Point right, double eps) {
	double Fy, Fz;
	Point y, z;
	Point length = right - left;
	Point lessEps = length / length.norm() * eps / 10;
	int funcCalls = 0;

	do {
		y = (left + right - lessEps) / 2;
		Fy = function(y);
		funcCalls++;

		z = (left + right + lessEps) / 2;
		Fz = function(z);
		funcCalls++;

		if (Fy < Fz) {
			right = z;
		}
		else {
			left = y;
		}
		length = right - left;
	} while (length.norm() > eps);

	Point ans = (right + left) / 2.;
	return ans;
}

int main() {
	Point point(0,0);
	Point grad;
	int iter = 0;

	cout << "Метод золотого сечения:\n\n";

	for (double eps = 0.1; eps > 0.00005; eps /= 10) {
		cout << "Точность:" << eps << endl;
		do {
			if (iter++ > 100) {
				cout << "Метод не смог дать ответа с данной точность за 1000 итераций.\n\n";
				return 0;
			}
			grad = gradient(point);
			point = GoldenRatio(point, point - grad, eps / 10);
		} while (grad.norm() > eps);

		cout << "Полученный ответ: (" << point.x << ", " << point.y << ")\n";
		cout << "Значение функции в данной точке: " << function(point) << "\n\n";

		iter = 0;
		point.x = 0;
		point.y = 0;
	}

	cout << "Метод дихотомии:\n\n";

	for (double eps = 0.1; eps > 0.00005; eps /= 10) {
		cout << "Точность:" << eps << endl;
		do {
			if (iter++ > 1000) {
				cout << "Метод не смог дать ответа с данной точность за 1000 итераций.";
				return 0;
			}
			grad = gradient(point);
			point = Dichotomy(point, point - grad, eps / 10);
		} while (grad.norm() > eps);

		cout << "Полученный ответ: (" << point.x << ", " << point.y << ")\n";
		cout << "Значение функции в данной точке: " << function(point) << "\n\n";

		iter = 0;
		point.x = 0;
		point.y = 0;
	}

	return 0;
}