#include <iostream>
#include "LPProblem.h"
#include "Point.h"
#include "math.h"
#include <vector>

# define M_PI 3.14159265358979323846

using namespace std;

double function(vector<double> x) {
	return x[0] * x[0] + 2 * x[1] * x[1] + sin(5 * x[0] + 3 * x[1]) + 3 * x[0] + 2 * x[1];
}

vector<double> gradient(vector<double> x) {
	vector<double> grad = { 2 * x[0] + 5 * cos(5 * x[0] + 3 * x[1]) + 3,
							4 * x[1] + 3 * cos(5 * x[0] + 3 * x[1]) + 2,
						   -1 };
	return grad;
}

bool is_in_omega(vector<double> x) { ///////////////
	return (6 * x[0] * x[0] + 9 * x[0] - 3 * x[1] <= 1) && (5 * x[0] + 3 * x[1] <= -2 * M_PI) && (x[0] - 5 * x[1] <= 3)
		&& (function(x) <= x[2]);
}

double vector_norm(vector<double> x) {
	return sqrt(x[0] * x[0] + x[1] * x[1]);
}

double phi_function(vector<double> x) {
	double max_value = 6 * x[0] * x[0] + 9 * x[0] - 3 * x[1] - 1;
	double temp = 5 * x[0] + 3 * x[1] + 2 * M_PI;
	if (temp > max_value)
		max_value = temp;
	temp = x[0] - 5 * x[1] - 3;
	if (temp > max_value)
		max_value = temp;
	temp = function(x) - x[2];
	if (temp > max_value)
		max_value = temp;
	return max_value;
}

vector<double> gradient_phi(vector<double> x) {
	int index = 0;
	double max_value = 6 * x[0] * x[0] + 9 * x[0] - 3 * x[1] - 1;
	double temp = 5 * x[0] + 3 * x[1] + 2 * M_PI;
	if (temp > max_value) {
		index = 1;
		max_value = temp;
	}
	temp = x[0] - 5 * x[1] - 3;
	if (temp > max_value) {
		index = 2;
		max_value = temp;
	}
	if (function(x) - x[2] > max_value)
		index = 3;

	vector<double> ans;
	if (index == 0) ans = { 12 * x[0] + 9, -3, 0 };
	if (index == 1) ans = { 5, 3, 0 };
	if (index == 2) ans = { -1, 2, 0 };
	if (index == 3) ans = gradient(x);

	return ans;
}

void cutting_hyperplane() {
	vector<double> fun_coeffs = { 0, 0, -1 };
	vector<vector<double>> res_matrix = { {5, 3, 0}, {1, -5, 0}, {-5, -3, 0}, {-1, 2, 0}, {3, 2, -1} };
	vector<double> res_vector = { -2 * M_PI, 3, 3 * M_PI, 4, 0 };
	vector<int> signs = { 1, 1, 1, 1, 1 };
	vector<int> res_signs = { 0, 0, 0 };
	double eps = 0.1;

	LPProblem* lp = new LPProblem(fun_coeffs, res_matrix, res_vector, signs, res_signs, 1);
	LPProblem* dual_lp = lp->dual();

	dual_lp->to_canonical_form();

	vector<double> x = dual_lp->simplex();
	vector<double> point = lp->direct_by_dual(x);
	vector<double> prev_point = point;
	vector<double> diff;
	int iter = 1;

	vector<double> a = gradient_phi(point);
	double b = (-1) * phi_function(point) + a[0] * point[0] + a[1] * point[1] + a[2] * point[2];
	res_matrix.push_back(a);
	res_vector.push_back(b);
	signs.push_back(1);

	while (eps > 0.00005) {
		lp = new LPProblem(fun_coeffs, res_matrix, res_vector, signs, res_signs, 1);

		dual_lp = lp->dual();
		dual_lp->to_canonical_form();

		x.push_back(0);
		x = dual_lp->simplex_method2(x);
		point = lp->direct_by_dual(x);
		diff = { point[0] - prev_point[0], point[1] - prev_point[1] };
		if (vector_norm(diff) < eps) {
			cout << "Точность: " << eps << endl;
			cout << "Точка: (" << point[0] << ", " << point[1] << "). x3 = " << point[2] << "\n";
			cout << "Норма градиента функции в точке: " << vector_norm(gradient(point)) << "\n";
			cout << "Значение функции в точке: " << function(point) << "\n";
			cout << "Количество итераций: " << iter << "\n\n";
			eps /= 10;
		}
		prev_point = point;
		a = gradient_phi(point);
		b = (-1) * phi_function(point) + a[0] * point[0] + a[1] * point[1] + a[2] * point[2];
		res_matrix.push_back(a);
		res_vector.push_back(b);
		signs.push_back(1);
		iter++;
	}
}

/* 
	Функция: 
		f(x) = x1^2 + 2x2^2 + sin(5x1 + 3x2) + 3x1 + 2x2
	Нелинейные ограничения:
		6x1^2 + 9x1 - 3x2 - 1 <= 0
		x1^2 + 2x2^2 + sin(5x1 + 3x2) + 3x1 + 2x2 - x3 <= 0
	Линейные ограничения:
		5x1 + 3x2 + 2pi <= 0
		x1 - 5x2 - 3 <= 0
	Дополнетельные линейные ограничения:
		-5x1 - 3x2 - 3pi <= 0
		-x1 + 2x2 - 4 <= 0
		3x1 + 2x2 - x3 - 1 <= 0
*/

int main() {
	cutting_hyperplane();
}