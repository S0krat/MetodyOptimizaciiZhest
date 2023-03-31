#include <iostream>
#include <vector>
#include "LPProblem.h"

using namespace std;

int main() {

	vector<double> f_coeffs = { 2, -1, -9, 1, -2 };

	vector<vector<double>> r_matrix = { {2,  -2,  5, -3,  1},
										{1, -1,  1,  2,  1},
										{2,  3, -3,  4,  2},
										{1, -1,  2,  1,  1},
										{3,  2,  1, -1,  3}, 
										{2, -2,  2, -2, -1},
										{4,  2, -3,  1, -3} };
	vector<double> r_vector = { 40, 30, 50, 25, 30, 45, 55 };

	vector<int> sign = { 0, 0, 0, 2, 2, 1, 1 };

	vector<int> r_signs = { 0, 0, 1, 1, 1};

	cout << "Прямая задача:\n\n";
	LPProblem* lp = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 1);
	lp->to_canonical_form();
	vector<double> x = lp->simplex();
	cout << "Ответ: " << x[0] - x[1] << " " << x[2] - x[3] << " " << x[4] << " " << x[5] << " " << x[6] << "\n";
	cout << "Целевая функция в точке: " << lp->obj_fun(x) << "\n\n";
	
	cout << "Двойственная задача:\n\n";
	LPProblem* lp2 = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 1);
	LPProblem* dual_lp = lp2->dual();
	dual_lp->to_canonical_form();
	x = dual_lp->extreme_point_enum();
	cout << "Ответ: " << x[0] - x[1] << " " << x[2] - x[3] << " " << x[4] - x[5] << " " << x[6] << " " << x[7] << " " << x[8] << " " << x[9] << "\n";
	cout << "Целевая функция в точке: " << dual_lp->obj_fun(x) << "\n\n";

	lp = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 1);
	vector<double> dual_sol = { x[0] - x[1], x[2] - x[3], x[4] - x[5], x[6], x[7], x[8], x[9] };
	x = lp->direct_by_dual(dual_sol);
	for (auto& cf : x) {
		cout << cf << " ";
	}
}