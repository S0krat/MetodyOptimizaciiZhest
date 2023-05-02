#include "TransportTask.h"

using namespace std;

int main() {
	vector<vector<int>> c = { {7, 1, 4, 14, 10},
							  {3, 7, 4,  9, 13},
							  {1, 6, 7,  8,  8},
							  {3, 2, 5,  6,  7} };
	vector<int> a = { 22, 17, 11, 12 };
	vector<int> b = { 16, 3, 11, 12, 20 };

	TransportTask tt = TransportTask(c, a, b);
	vector<vector<int>> sol = tt.northwest_solution();

	/*for (auto& row : sol) {
		for (int& elem : row) {
			cout << elem << " ";
		}
		cout << endl;
	}*/

	tt.potential_method(sol);
}