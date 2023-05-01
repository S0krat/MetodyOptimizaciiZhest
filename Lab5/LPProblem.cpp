#include "LPProblem.h"

void LPProblem::vector_print(vector<double>& v) {
	for (auto& elem : v) {
		cout << elem << " ";
	}
	cout << "\n";
}

void LPProblem::matrix_print(vector<vector<double>>& a) {
	for (int i = 0; i < a.size(); i++) {
		for (auto& elem : a[i]) {
			if (elem >= 0) cout << " ";
			cout << elem << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}

void LPProblem::print() {
	for (int i = 0; i < var_num; i++) {
		if (i != 0) cout << " + ";
		cout << fun_coeffs[i] << "x_" << i + 1;
	}
	cout << " --> " << (cond ? "max\n\n" : "min\n\n");

	for (int i = 0; i < res_num; i++) {
		for (int j = 0; j < var_num; j++) {
			if (j != 0) cout << " + ";
			cout << res_matrix[i][j] << "x_" << j + 1;
		}

		if (signs[i] == 0) cout << " = ";
		else if (signs[i] == 1) cout << " <= ";
		else if (signs[i] == 2) cout << " >= ";

		cout << res_vector[i] << endl;
	}

	for (int i = 0; i < var_num; i++) {
		if (res_signs[i] == 1) cout << "x_" << i + 1 << " >= 0; ";
		if (res_signs[i] == 2) cout << "x_" << i + 1 << " <= 0; ";
	}
	cout << "\n\n";
}

void LPProblem::to_canonical_form() {
	if (cond == 1) {
		for (int i = 0; i < var_num; i++)
			fun_coeffs[i] *= -1;
		cond = 0;
		condChange = true;
	}

	for (int i = 0; i < res_num; i++) {
		if (res_vector[i] < 0) {
			res_vector[i] *= -1;
			for (int j = 0; j < var_num; j++) {
				res_matrix[i][j] *= -1;
			}
			if (signs[i] == 1) signs[i] = 2;
			else if (signs[i] == 2) signs[i] = 1;
		}
	}

	for (int i = 0; i < var_num; i++) {
		if (res_signs[i] == 2) {
			res_signs[i] = 1;
			fun_coeffs[i] *= -1;
			for (int j = 0; j < res_num; j++) {
				res_matrix[j][i] *= -1;
			}
		}
		else if (res_signs[i] == 0) {
			res_signs[i] = 1;
			for (int j = 0; j < res_num; j++) {
				res_matrix[j].push_back(res_matrix[j][var_num - 1]);
				for (int k = var_num - 1; k > i + 1; k--) {
					res_matrix[j][k] = res_matrix[j][k - 1];
				}
			}
			res_signs.push_back(res_signs[var_num - 1]);
			for (int j = var_num - 1; j > i + 1; j--) {
				res_signs[j] = res_signs[j - 1];
			}
			res_signs[i + 1] = 1;
			res_signs[i] = 1;

			fun_coeffs.push_back(fun_coeffs[var_num - 1]);
			for (int j = var_num - 1; j > i + 1; j--) {
				fun_coeffs[j] = fun_coeffs[j - 1];
			}
			fun_coeffs[i + 1] = (-1) * fun_coeffs[i];

			var_num++;
			for (int j = 0; j < res_num; j++) {
				res_matrix[j][i + 1] = (-1) * res_matrix[j][i];
			}
		}

	}


	for (int i = 0; i < res_num; i++) {
		if (signs[i]) {
			fun_coeffs.push_back(0);
			for (int j = 0; j < res_num; j++) {
				if (i == j) res_matrix[j].push_back((signs[i] == 1 ? 1 : -1));
				else res_matrix[j].push_back(0);
			}
			signs[i] = 0;
			res_signs.push_back(1);
			var_num++;
		}
	}
}

vector<double> LPProblem::simplex_method2(vector<double>& sol) {
	if (sol.size() != var_num) exit(-2);
	const double zero = 0.000001;

	vector<int> Nk, Lk;
	for (int i = 0; i < var_num; i++) {
		if (sol[i] > zero) Nk.push_back(i);
		else Lk.push_back(i);
	}

	while (true) {
		bool isUnpl = true;
		vector<int> zeroCfs;
		if (Nk.size() > res_num) {
			cout << "Блииин...\n";
			exit(-3);
		}
		while (Nk.size() < res_num) {
			isUnpl = false;
			zeroCfs.push_back(Nk.size());
			Nk.push_back(Lk[Lk.size() - 1]);
			Lk.pop_back();
		}

		vector<vector<double>> ANk(res_num);
		vector<vector<double>> ALk(res_num);
		for (int i = 0; i < res_num; i++) {
			for (int& cf : Nk) {
				ANk[i].push_back(res_matrix[i][cf]);
			}
			for (int& cf : Lk) {
				ALk[i].push_back(res_matrix[i][cf]);
			}
		}

		vector<vector<double>> B = inverse_matrix(ANk);
		if (B.size() == 1) {
			Lk.push_back(Nk[Nk.size() - 1]);
			Nk.pop_back();
			Nk.push_back(Lk[0]);
			Lk.erase(Lk.begin());
			continue;
		}

		vector<double> cLk;
		vector<double> cNk;
		for (int& cf : Lk) {
			cLk.push_back(fun_coeffs[cf]);
		}
		for (int& cf : Nk) {
			cNk.push_back(fun_coeffs[cf]);
		}

		vector<double> cNkBNk(res_num);
		for (int i = 0; i < res_num; i++) {
			for (int j = 0; j < res_num; j++) {
				cNkBNk[i] += cNk[j] * B[j][i];
			}
		}

		vector<double> cNkBNkALk(var_num - res_num);
		for (int i = 0; i < var_num - res_num; i++) {
			for (int j = 0; j < res_num; j++) {
				cNkBNkALk[i] += cNkBNk[j] * ALk[j][i];
			}
		}

		vector<double> dLk(var_num - res_num);
		for (int i = 0; i < var_num - res_num; i++) {
			dLk[i] = cLk[i] - cNkBNkALk[i];
		}

		bool isOpt = true;
		int jk = 0;
		for (int i = 0; i < var_num - res_num; i++) {
			if (dLk[i] <= -zero) {
				isOpt = false;
				jk = i;
			}
		}
		if (isOpt) {
			for (int i = 0; i < var_num; i++) {
				if (sol[i] < zero) sol[i] = 0;
			}
			return sol;
		}

		vector<double> uk(res_num);
		for (int i = 0; i < res_num; i++) {
			for (int j = 0; j < res_num; j++) {
				uk[i] += B[i][j] * ALk[j][jk];
			}
		}
		
		bool isInf = true;
		vector<int> ik;
		for (int i = 0; i < res_num; i++) {
			if (uk[i] > 0) {
				isInf = false;
				ik.push_back(i);
			}
		}
		if (isInf) {
			cout << "Множество допустимых решений неограниченно.";
			exit(-4);
		}

		if (isUnpl) {
			double theta = sol[Nk[ik[0]]] / uk[ik[0]];
			for (int i = 1; i < ik.size(); i++) {
				double temp = sol[Nk[ik[i]]] / uk[ik[i]];
				if (temp < theta && temp > zero) theta = temp;
			}

			for (int i = 0; i < res_num; i++) {
				sol[Nk[i]] -= theta * uk[i];
			}
			sol[Lk[jk]] += theta;

			Lk.clear();
			Nk.clear();
			for (int i = 0; i < var_num; i++) {
				if (sol[i] > zero) Nk.push_back(i);
				else Lk.push_back(i);
			}
		}
		else {
			bool isNeg = true;
			for (int& cf : zeroCfs) {
				if (uk[cf] > zero) isNeg = false;
			}
			if (isNeg) {
				double theta = INT_MAX;
				for (int i = 0; i < ik.size(); i++) {
					double temp = sol[Nk[ik[i]]] / uk[ik[i]];
					if (temp < theta && temp > zero) theta = temp;
				}

				for (int i = 0; i < res_num; i++) {
					sol[Nk[i]] -= theta * uk[i];
				}
				sol[Lk[jk]] += theta;
				Lk.clear();
				Nk.clear();
				for (int i = 0; i < var_num; i++) {
					if (sol[i] > zero) Nk.push_back(i);
					else Lk.push_back(i);
				}
			}
			else {
				Lk.push_back(Nk[Nk.size() - 1]);
				Nk.pop_back();
				Nk.push_back(Lk[0]);
				Lk.erase(Lk.begin());
			}
		}
	}
}

vector<double> LPProblem::simplex() {
	vector<vector<double>> A = res_matrix;
	for (int i = 0; i < res_num; i++) {
		for (int j = 0; j < res_num; j++) {
			if (i == j) A[i].push_back(1);
			else A[i].push_back(0);
		}
	}
	vector<double> c;
	
	for (int i = 0; i < var_num; i++) {
		c.push_back(0);
	}
	for (int i = 0; i < res_num; i++) {
		c.push_back(1);
	}
	vector<double> r_vec = res_vector;
	vector<int> s = { 1,1 };
	vector<int> r_s = { 1,1 };

	LPProblem* new_lp = new LPProblem(c, A, r_vec, s, r_s, 0);
	vector<double> acc_sol;
	for (int i = 0; i < var_num; i++) {
		acc_sol.push_back(0);
	}
	for (int i = 0; i < res_num; i++) {
		acc_sol.push_back(res_vector[i]);
	}
	vector<double> solution = new_lp->extreme_point_enum();

	for (int i = 0; i < res_num; i++) {
		solution.pop_back();
	}
	//cout << "Найдено допустимое решение: ";
	//this->vector_print(solution);
	vector<double> sol = this->simplex_method2(solution);
	return sol;
}

LPProblem* LPProblem::dual() {
	int new_var_num = res_num;
	int new_res_num = var_num;
	vector<double> new_fun_coeffs = res_vector;
	vector<double> new_res_vector = fun_coeffs;
	vector<int> new_signs;
	for (auto& sign : res_signs) {
		if (sign == 0) new_signs.push_back(0);
		else if (sign == 1) new_signs.push_back(2);
		else if (sign == 2) new_signs.push_back(1);
	}
	vector<int> new_res_signs = signs;
	vector<vector<double>> new_res_matrix(new_res_num);
	for (int i = 0; i < new_res_num; i++) {
		for (int j = 0; j < new_var_num; j++) {
			new_res_matrix[i].push_back(res_matrix[j][i]);
		}
	}
	int new_cond = (cond ? 0 : 1);
	return new LPProblem(new_fun_coeffs, new_res_matrix, new_res_vector, new_signs, new_res_signs, new_cond);
}

vector<double> LPProblem::extreme_point_enum() {
	string bitmask(res_num, 1);
	bitmask.resize(var_num, 0);
	vector<vector<double>> matrix(res_num);
	for (int j = 0; j < res_num; j++) {
		for (int i = 0; i < res_num; i++) {
			matrix[j].push_back(0);
		}
	}
	vector<double> solution(res_num);
	vector<double> best_solution(res_num);
	vector<int> vars(res_num);
	vector<int> vars_best(res_num);
	double min_sum = INT_MAX;

	do {
		int coef = 0;
		for (int i = 0; i < var_num; ++i) {
			if (bitmask[i]) {
				for (int j = 0; j < res_num; j++) {
					matrix[j][coef] = res_matrix[j][i];
				}
				vars[coef] = i;
				coef++;
			}
			
		}

		solution = gauss(matrix, res_vector, res_num);

		double sum = 0;
		bool isAppr = false;
		for (int i = 0; i < res_num; i++) {
			if (solution[i] != 0) isAppr = true;
		}
		for (int i = 0; i < res_num; i++) {
			if (solution[i] < 0) isAppr = false;
			sum += solution[i] * fun_coeffs[vars[i]];
		}


		//if (isAppr) {
		//	for (int i = 0; i < res_num; i++) {
		//		cout << "x_" << vars[i] + 1 << " = " << solution[i] << "; ";
		//	}
		//	cout << "\nFunction = " << sum;
		//	cout << "\n";
		//}

		if (isAppr && sum < min_sum) {
			min_sum = sum;
			best_solution = solution;
			vars_best = vars;
			//cout << "Есть допустимое решение!\n";
		}
	} while (prev_permutation(bitmask.begin(), bitmask.end()));

	vector<double> sol(var_num);
	for (int i = 0; i < res_num; i++) {
		sol[vars_best[i]] = best_solution[i];
	}

	return sol;
}

vector<double> LPProblem::direct_by_dual(vector<double>& sol) {
	if (sol.size() != res_num)
		exit(-1);
	double eps = 0.00001;

	vector<int> nonzero;
	for (int i = 0; i < var_num; i++) {
		double sum = 0;
		for (int j = 0; j < res_num; j++) {
			sum += sol[j] * res_matrix[j][i];
		}
		if (abs(sum - fun_coeffs[i]) < eps) {
			nonzero.push_back(i);
		}
	}

	vector<vector<double>> matrix(nonzero.size());
	vector<double> vect(nonzero.size());
	int coef = 0;
	for (int i = 0; i < res_num; i++) {
		if (abs(sol[i]) > eps) {
			for (auto& cf : nonzero) {
				matrix[coef].push_back(res_matrix[i][cf]);
			}
			vect[coef] = res_vector[i];
			coef++;
		}
	}
	
	vector<double> gaus = gauss(matrix, vect, nonzero.size());
	vector<double> ans(var_num);
	coef = 0;
	for (auto& cf : nonzero) {
		ans[cf] = gaus[coef++];
	}
	
	return ans;
}

vector<double> LPProblem::gauss(vector<vector<double>>& A, vector<double>& Y, int n) {
	vector<double> x(n);
	vector<vector<double>> a = A;
	vector<double> y = Y;
	double max;
	int k, index;
	const double eps = 0.00001;  // точность
	k = 0;
	while (k < n) {
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++) {
			if (abs(a[i][k]) > max) {
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps) {
			// нет ненулевых диагональных элементов
			//cout << "Решение получить невозможно из-за нулевого столбца ";
			//cout << index << " матрицы A" << endl;
			return x;
		}
		for (int j = 0; j < n; j++) {
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++) {
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}

vector<vector<double>> LPProblem::inverse_matrix(vector<vector<double>>& A) {
	int n = A.size();
	vector<vector<double>> E(n);
	vector<vector<double>> a = A;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			E[i].push_back(((i == j) ? 1 : 0));
		}
	}
	double max;
	int k, index;
	const double eps = 0.00001;  // точность
	k = 0;
	while (k < n) {
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++) {
			if (abs(a[i][k]) > max) {
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps) {
			// нет ненулевых диагональных элементов
			//cout << "Обратная матрица не существует ";
			//cout << index << endl;
			return { {0} };
		}
		for (int j = 0; j < n; j++) {
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
			temp = E[k][j];
			E[k][j] = E[index][j];
			E[index][j] = temp;
		}
		// Нормализация уравнений
		for (int i = k; i < n; i++) {
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++) {
				a[i][j] = a[i][j] / temp;
				E[i][j] = E[i][j] / temp;
			}
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++) {
				a[i][j] = a[i][j] - a[k][j];
				E[i][j] = E[i][j] - E[k][j];
			}		
		}
		k++;
	}
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			for (int k = 0; k < n; k++) {
				E[j][k] -= E[i][k] * a[j][i];
			}
		}
	}
	return E;
}