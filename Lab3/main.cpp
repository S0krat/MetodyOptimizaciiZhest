#include <iostream>
#include <math.h>

using namespace std;

double function(double x) {
	return (x - 1) / (x * x);
}

void TrialPoint(double left, double right, double eps) {
	double mid, y, z, Fmid, Fy, Fz;
	double length = right - left;
	int funcCalls = 0;

	mid = (left + right) / 2.;
	Fmid = function(mid);
	funcCalls++;

	do {
		y = left + length / 4;
		Fy = function(y);
		funcCalls++;

		z = right - length / 4;
		Fz = function(z);
		funcCalls++;

		if (Fy > Fmid) {
			right = mid;
			mid = y;
			Fmid = Fy;
			length /= 2.;
			continue;
		}
		if (Fz > Fmid) {
			left = mid;
			mid = z;
			Fmid = Fz;
		}
		else {
			left = y;
			right = z;
		}
		length /= 2.;
	} while (length > eps);
	
	cout << "Метод пробных точек:\nОтвет: " << (right + left) / 2. << " +- " << eps / 2 << endl;
	cout << "Количество обращений к функции: " << funcCalls << "\n\n";
}

void Dichotomy(double left, double right, double eps) {
	double y, z, Fy, Fz;
	double length = right - left;
	double lessEps = eps / 10;
	int funcCalls = 0;

	do {
		y = (left + right - lessEps) / 2;
		Fy = function(y);
		funcCalls++;

		z = (left + right + lessEps) / 2;
		Fz = function(z);
		funcCalls++;

		if (Fy > Fz) {
			right = z;
		}
		else {
			left = y;
		}
		length = right - left;
	} while (length > eps);

	cout << "Метод дихотомии:\nОтвет: " << (right + left) / 2. << " +- " << eps / 2 << endl;
	cout << "Количество обращений к функции: " << funcCalls << "\n\n";
}

void GoldenRatio(double left, double right, double eps) {
	const double phi = (3 - sqrt(5)) / 2;
	double y, z, Fy, Fz;
	int funcCalls = 0;

	y = left + phi * (right - left);
	Fy = function(y);
	funcCalls++;

	z = right + left - y;
	Fz = function(z);
	funcCalls++;

	while (true) {
		if (Fy > Fz) {
			right = z;
			if (right - left < eps) break;
			z = y;
			Fz = Fy;
			y = right + left - y;
			Fy = function(y);
			funcCalls++;
		} else {
			left = y;
			if (right - left < eps) break;
			y = z;
			Fy = Fz;
			z = left + right - z;
			Fz = function(z);
			funcCalls++;
		}
	}

	cout << "Метод золотого сечения:\nОтвет: " << (right + left) / 2. << " +- " << eps / 2 << endl;
	cout << "Количество обращений к функции: " << funcCalls << "\n\n";
}

int main() {
	cout.setf(ios_base::fixed);
	cout.precision(7);
	const double LEFT = 0.1;
	const double RIGHT = 2.9;
	const double epsilon = 0.000001;
	TrialPoint(LEFT, RIGHT, epsilon);
	Dichotomy(LEFT, RIGHT, epsilon);
	GoldenRatio(LEFT, RIGHT, epsilon);
}