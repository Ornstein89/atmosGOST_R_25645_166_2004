/* Test and example of use of the atmosGOST_R_25645_166_2004() function */

#include <iostream>
#include <chrono>

#include "atmosGOST_R_25645_166_2004.h"

double F0_arr_test[7] = { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };

double rand_f(double a, double b) {
	return ((double)rand() / RAND_MAX) * (b - a) + a;
}

int main()
{
	// TODO тесты по таблице

	/*
	// тест быстродействия алгоритма поиска ближайшего
	srand(time(NULL));
	// алгоритм 1 на основе прохода по таблице
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; i++)
	{
		double F81 = rand_f(65.0, 270.0);
		int n_col = 0;
		double diff1 = abs(F0_arr_test[n_col] - F81);
		double diff2 = diff1;
		while (n_col < 6)
		{
			n_col++;
			diff2 = abs(F0_arr_test[n_col] - F81);
			if (diff1 <= diff2) {
				--n_col;
				break;
			}
			diff1 = diff2;
		}
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Time 1: " << elapsed.count() << std::endl;

	// алгоритм 2 на основе деления
	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; i++)
	{
		double F81 = rand_f(65.0, 270.0);
		double F0 = 0.0;
		int n_col = 0.0;
		if (F81 <= 200.0 && F81 >= 75.0) {
			double tmp = round(F81 / 25.0);
			n_col = tmp - 3;
			F0 = tmp * 25.0;
		}
		else if (F81 <= 250.0 && F81 >= 200.0)
		{
			double tmp = round(F81 / 50.0);
			n_col = tmp + 1;
			F0 = tmp * 50.0;
		}
		else if (F81 < 75.0) {
			F0 = 75.0;
			n_col = 0;
		}
		else {
			F0 = round(F81 / 25) * 25.0;
			n_col = 6;
		}
	}
	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	std::cout << "Time 2: " << elapsed.count() << std::endl;
	*/

	std::system("pause");
	return 0;
}