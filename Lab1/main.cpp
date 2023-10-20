#include <iostream>
#include <vector>
#include "LPProblem.h"

using namespace std;

int main() {
	setlocale(LC_ALL, "Russian");
	vector<double> f_coeffs = { 3, 4, -2 };

	vector<vector<double>> r_matrix = { {-1,  1,  1},
										{-2, 3,  -1},
										{-1,  -2, 4} };
	vector<double> r_vector = { 10, 5, 20 };

	vector<int> sign = { 0, 2, 1 };

	vector<int> r_signs = { 1, 1, 0};

	cout << "Прямая задача:\n\n";
	LPProblem* lp = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 0);
	lp->print();
	lp->to_canonical_form();
	vector<double> x = lp->simplex();
	cout << "Ответ: " << x[0] << " " << x[1] << " " << x[2]-x[3] << "\n";
	cout << "Целевая функция в точке: " << lp->obj_fun(x) << "\n\n";
	
	cout << "Двойственная задача:\n\n";
	LPProblem* lp2 = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 0);
	LPProblem* dual_lp = lp2->dual();
	dual_lp->print();
	dual_lp->to_canonical_form();
	dual_lp->print();
	x = dual_lp->simplex();
	cout << "Ответ: " << x[0] << " " << x[1] << " " << x[2] - x[3] << "\n";
	cout << "Целевая функция в точке: " << dual_lp->obj_fun(x) << "\n\n";
}