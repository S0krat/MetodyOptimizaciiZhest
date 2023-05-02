#pragma once
#include <vector>
#include <iostream>

using namespace std;

class TransportTask {
public:
	TransportTask(vector<vector<int>>& c, vector<int>& a, vector<int>& b) {
		C = c;
		A = a;
		M = a.size();
		B = b;
		N = b.size();
	}
	vector<vector<int>> northwest_solution() {
		vector<int> A_temp = A;
		vector<int> B_temp = B;
		vector<vector<int>> sol(M);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				sol[i].push_back(0);
			}
		}
 		int i = 0; 
		int j = 0;
		while (i != M && j != N) {
			if (A_temp[i] < B_temp[j]) {
				sol[i][j] = A_temp[i];
				B_temp[j] -= A_temp[i];
				A_temp[i] = 0;
				i++;
			}
			else {
				sol[i][j] = B_temp[j];
				A_temp[i] -= B_temp[j];
				B_temp[j] = 0;
				j++;
			}
		}
		return sol;
	}
	vector<vector<int>> potential_method(vector<vector<int>>& sol) {
		vector<int> u(M);
		vector<int> v(N);
		v[0] = C[0][0];
		int i = 0;
		int j = 0;
		while (i < M - 1 || j < N - 1) {
			if (i < M - 1 && sol[i + 1][j] > 0) {
				i++;
				u[i] = v[j] - C[i][j];
			}
			else {
				j++;
				v[j] = u[i] + C[i][j];
			}
		}
		//cout << "u === ";
		//for (auto& row : u) {
		//	cout << row << " ";
		//}
		//cout << "\nv === ";
		//for (auto& row : v) {
		//	cout << row << " ";
		//}
		//cout << "\n\n";

		vector<vector<int>> alpha(M);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				alpha[i].push_back(v[j] - u[i]);
			}
		}

		while (true) {
			bool isOpt = true;
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < N; j++) {
					if (alpha[i][j] > C[i][j]) isOpt = false;
				}
			}
			if (isOpt) {
				return sol;
			}
		}
		//for (auto& row : alpha) {
		//	for (int& elem : row) {
		//		cout << elem << " ";
		//	}
		//	cout << endl;
		//}
		return sol;
	}
private:
	int M;
	int N;
	vector<vector<int>> C;
	vector<int> A;
	vector<int> B;
};

