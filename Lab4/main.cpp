#include "Point.h"
#include <iostream>

using namespace std;

double function(Point p) {
	return p.x * p.x + 2 * p.y * p.y + sin(5 * p.x + 3 * p.y) + 3 * p.x + 2 * p.y;
}

Point gradient(Point p) {
	return Point(2 * p.x + 5 * cos(5 * p.x + 3 * p.y) + 3, 4 * p.y + 3 * cos(5 * p.x + 3 * p.y) + 2);
}

Matrix Hesse(Point p) {
	return Matrix(Point(2 - 25 * sin(5 * p.x + 3 * p.y), -15 * sin(5 * p.x + 3 * p.y)), 
				  Point(-15 * sin(5 * p.x + 3 * p.y), 4 - 9 * sin(5 * p.x + 3 * p.y)));
}

Point GoldenRatio(Point left, Point right, double eps, int* funcCalls) {
	const double phi = (3 - sqrt(5)) / 2;
	Point y, z;
	double Fy, Fz;

	y = right - left;
	y = y * phi;
	y = y + left;

	Fy = function(y);
	(*funcCalls)++;

	z = right + (left - y);
	Fz = function(z);
	(*funcCalls)++;

	while (true) {
		if (Fy < Fz) {
			right = z;
			if (right.distTo(left) < eps) break;
			z = y;
			Fz = Fy;
			y = right + (left - y);
			Fy = function(y);
			(*funcCalls)++;
		}
		else {
			left = y;
			if (right.distTo(left) < eps) break;
			y = z;
			Fy = Fz;
			z = left + (right - z);
			Fz = function(z);
			(*funcCalls)++;
		}
	}

	Point ans = (right + left) / 2.;
	return ans;
}

Point Dichotomy(Point left, Point right, double eps, int* funcCalls) {
	double Fy, Fz;
	Point y, z;
	Point length = right - left;
	Point lessEps = length / length.norm() * eps / 10;

	do {
		y = (left + right - lessEps) / 2;
		Fy = function(y);
		(*funcCalls)++;

		z = (left + right + lessEps) / 2;
		Fz = function(z);
		(*funcCalls)++;

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

Point getPoint(Point start, Point grad) {
	Point sumPoint = start - grad;
	double sumPointValue = 5 * sumPoint.x + 3 * sumPoint.y;
	if (sumPointValue > -6.283) {
		double startPointValue = 5 * start.x + 3 * start.y;
		grad = grad * ((-6.283 - startPointValue) / (sumPointValue - startPointValue));
		return start - grad;
	}
	else if (sumPointValue < -9.4248) {
		double startPointValue = 5 * start.x + 3 * start.y;
		grad = grad * ((startPointValue + 9.4248) / (startPointValue - sumPointValue));
		return start - grad;
	}
	return sumPoint;
}

void function1() {
	Point x_star = Point(-1.3186799325, -0.4456972731);

	const double start_p1 = -1;
	const double start_p2 = -1;
	Point point(start_p1, start_p2);
	Point tempPoint;
	Point grad;
	Point prevPoint;
	int iter = 0;
	int funcCalls = 0;
	double coef = 2 * (1 + 2 / 27);

		do {
			if (iter++ > 1000) {
				cout << "Метод не смог дать ответа с данной точность за 1000 итераций.\n\n";
				exit(0);
			}
			grad = gradient(point);

			// cout << "Point: " << point.x << " " << point.y << endl;
			tempPoint = getPoint(point, grad);
			prevPoint = point;
			point = GoldenRatio(point, tempPoint, 0.001, &funcCalls);
			grad = gradient(point);
			cout << grad.norm() * grad.norm() << "   ";
			cout << coef * (function(point) - function(x_star)) << "   ";
			grad = gradient(prevPoint);
			cout << (prevPoint.x - point.x) / grad.x << "\n";
		} while (grad.norm() > 0.01);

		cout << 2. / 27. * 0.99;
}

void function2() {
	Point x_star = Point(-1.31810253653345932356, -0.44543076096003786368);

	const double start_p1 = -1;
	const double start_p2 = -1;
	Point point(start_p1, start_p2);
	Point tempPoint;
	Point grad;
	Point temp1, temp2;
	int iter = 0;
	int funcCalls = 0;
	Matrix hesse;
	Point gradhesse;
	double value, norm1, norm2;

		do {
			if (iter++ > 100) {
				cout << "Метод не смог дать ответа с данной точность за 1000 итераций.";
				exit(0);
			}
			grad = gradient(point);
			value = function(point);
			funcCalls++;
			hesse = Hesse(point);
			hesse.inverse();
			gradhesse = (hesse * grad) * (-1);
			while (function(point + gradhesse) - value > (grad * gradhesse) / 3) {
				gradhesse = gradhesse * 0.5;
				funcCalls++;
			}
			temp1 = point - x_star;
			norm1 = temp1.norm();
			//cout << "||x_k - x_*|| = " << norm1 << endl;
			point = point + gradhesse;
			//cout << "Точка: " << point.x << " " << point.y << endl;
			temp2 = point - x_star;
			norm2 = temp2.norm();
			//cout << "||x_(k+1) - x_*|| = " << norm2 << endl;
			cout << norm2 / (norm1 * norm1) << "\n\n";
		} while (grad.norm() > 0.0001);
}

int main() {
	//cout.setf(ios_base::fixed);
	//cout.precision(20);
	//function2();
	//return 0;
	const double start_p1 = -1;
	const double start_p2 = -1;
	Point point(start_p1, start_p2);
	Point tempPoint;
	Point grad;
	int iter = 0;
	int funcCalls = 0;
	cout << "Методы первого порядка:\n";
	cout << "Метод золотого сечения:\n\n";

	for (double eps = 0.1; eps > 0.00005; eps /= 10) {
		cout << "Точность:" << eps << endl;
		do {
			if (iter++ > 1000) {
				cout << "Метод не смог дать ответа с данной точность за 1000 итераций.\n\n";
				return 0;
			}
			grad = gradient(point);
			// cout << "Point: " << point.x << " " << point.y << endl;
			tempPoint = getPoint(point, grad);
			point = GoldenRatio(point, tempPoint, eps / 5, &funcCalls);
		} while (grad.norm() > eps);

		cout << "Полученный ответ: (" << point.x << ", " << point.y << ")\n";
		cout << "Значение функции в данной точке: " << function(point) << "\n";
		cout << "Колчество вызовов функции: " << funcCalls << "\n";
		cout << "Количество вызовов градиента: " << iter << "\n\n";
		cout << "Норма градиента функции:" << grad.norm() << "\n";
		/*Matrix hesse = Hesse(point);
		cout << "Определитель гессиана в точке: " << hesse.det() << "\n\n";*/

		iter = 0;
		funcCalls = 0;
		point.x = start_p1;
		point.y = start_p2;
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
			tempPoint = getPoint(point, grad);
			point = Dichotomy(point, tempPoint, eps / 5, &funcCalls);
		} while (grad.norm() > eps);

		cout << "Полученный ответ: (" << point.x << ", " << point.y << ")\n";
		cout << "Значение функции в данной точке: " << function(point) << "\n";
		cout << "Колчество вызовов функции: " << funcCalls << "\n";
		cout << "Количество вызовов градиента: " << iter << "\n\n";
		cout << "Норма градиента функции:" << grad.norm() << "\n";
		/*Matrix hesse = Hesse(point);
		cout << "Определитель гессиана в точке: " << hesse.det() << "\n\n";*/

		iter = 0;
		funcCalls = 0;
		point.x = start_p1;
		point.y = start_p2;
	}

	cout << "Метод второго порядка:\n\n";
	Matrix hesse;
	Point gradhesse;
	double value;

	for (double eps = 0.1; eps > 0.00000005; eps /= 10) {
		cout << "Точность:" << eps << endl;
		do {
			if (iter++ > 100) {
				cout << "Метод не смог дать ответа с данной точность за 1000 итераций.";
				return 0;
			}
			grad = gradient(point);
			value = function(point);
			funcCalls++;
			hesse = Hesse(point);
			hesse.inverse();
			gradhesse = hesse * grad * (-1);
			while (function(point + gradhesse) - value > (grad * gradhesse) / 3) {
				gradhesse = gradhesse * 0.5;
				funcCalls++;
			}
			point = point + gradhesse;
		} while (grad.norm() > eps);

		cout << "Полученный ответ: (" << point.x << ", " << point.y << ")\n";
		cout << "Значение функции в данной точке: " << function(point) << "\n";
		cout << "Колчество вызовов функции: " << funcCalls << "\n";
		cout << "Количество вызовов градиента и гессиана: " << iter << "\n\n";
		/*cout << "Норма градиента функции:" << grad.norm() << "\n";
		hesse = Hesse(point);
		cout << "Определитель гессиана в точке: " << hesse.det() << "\n\n";*/

		iter = 0;
		funcCalls = 0;
		point.x = start_p1;
		point.y = start_p2;
	}

	return 0;
}