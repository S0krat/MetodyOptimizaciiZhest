#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>

using namespace std;

class LPProblem {
public:

	LPProblem(vector<double>& f_coeffs, vector<vector<double>>& r_matrix, vector<double>& r_vector, vector<int>& sign, vector<int>& r_signs, bool con) {
		var_num = f_coeffs.size();
		fun_coeffs = f_coeffs;
		res_num = r_vector.size();
		res_matrix = r_matrix;
		res_vector = r_vector;
		signs = sign;
		res_signs = r_signs;
		cond = con;
		condChange = false;
	}

	void print(); // Выводит задачу в консоль

	void to_canonical_form(); // Переписывает задачу в каноническую форму

	LPProblem* dual(); // Возвращает новую задачу - двойственную к данной

	vector<double> simplex_method2(vector<double>& sol);
	vector<double> simplex();

	vector<double> extreme_point_enum(); // Решает задачу методом перебора крайних точек

	vector<double> direct_by_dual(vector<double>& sol);

	void matrix_print(vector<vector<double>>& a);
	void vector_print(vector<double>& v);

	double obj_fun(vector<double> x) {
		double sum = 0;
		for (int i = 0; i < var_num; i++) {
			sum += x[i] * fun_coeffs[i];
		}
		return condChange ? -sum : sum;
	}
private:
	static vector<double> gauss(vector<vector<double>>& A, vector<double>& Y, int n);
	vector<vector<double>> inverse_matrix(vector<vector<double>>& A);
	// Количество переменных
	int var_num;
	// коэффициенты целевой функции
	vector<double> fun_coeffs;
	// количество ограничений
	int res_num;
	// матрица коэффициентов ограничений
	vector<vector<double>> res_matrix;
	// правый вектор-столбец ограничений
	vector<double> res_vector;
	// знаки в ограниениях (0 это '=', 1 это '<=', 2 это '>=')
	vector<int> signs;
	// ограничения на знак (0 - нет, 1 - '>= 0', 2 - '<= 0')
	vector<int> res_signs;
	// условие задачи (0 - минимизация цф, 1 - максимизация цф)
	bool cond;
	// менялось ли условие задачи
	bool condChange;
};