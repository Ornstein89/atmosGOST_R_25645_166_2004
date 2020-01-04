/* Test and example of use of the atmosGOST_R_25645_166_2004() function */

#include <iostream>
#include <string>
#include <chrono>

#include "atmosGOST_R_25645_166_2004.h"
#include "tables_4_9_GOST_R_25645_166_2004.h"

double F0_arr_test[7] = { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };

double rand_f(double a, double b) {
	return ((double)rand() / RAND_MAX) * (b - a) + a;
}

bool test_by_tables4_9(
	double const tol_rho,
	double const tol_coef,
	bool const verbose)
{
	double max_rho_err = 0.0;
	double max_coef_err = 0.0;
	for (int i = 0; i < 7; i++) // перечисление по шкале F
	{
		double max_rho_err_i = 0.0;
		double max_coef_err_i = 0.0;
		for (int j = 0; j < 70; j++)
		{
			double h = H0 + j * H_step;
			double X0[3] = { 0.0 };
			X0[0] = h;

			double rho = atmosGOST_R_25645_166_2004(
				h,
				F0_arr_test[i],
				1.0,
				F0_arr_test[i], // 
				0.0,
				X0,
				0.0,
				0.0,
				0.0,
				0.0);

			double err_rho0 = abs(test_glob_rho_night - table_rho_night[j][i]);
			max_rho_err_i = std::fmax(max_rho_err_i, err_rho0);

			double err_rhoK0 = abs(test_glob_K0 - table_K0[j][i]);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK0);
			double err_rhoK1 = abs(test_glob_K1 - table_K1[j][i]);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK1);
			double err_rhoK2 = abs(test_glob_K2 - table_K2[j][i]);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK2);
			double err_rhoK3 = abs(test_glob_K3 - table_K3[j][i]);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK3);
			double err_rhoK4_1 = abs(test_glob_K4_1 - table_K4[j][i]);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK4_1);
			
			//test_glob_K4_2 - ;

			if (verbose) {
				std::string rslt =
					((max_rho_err_i < tol_rho)&&(max_coef_err_i < tol_coef))
					? ("OK")
					: ("***FAILED");
				std::cout << rslt
					<< "\th = " << h
					<< ";\tF107 = " << F0_arr_test[i]
					<< ";\tmax_rho_error = " << max_rho_err_i
					<< ";\tmax_coef_error = " << max_coef_err_i
					<< std::endl;
			}
		}
		max_rho_err = std::fmax(max_rho_err, max_rho_err_i);
		max_coef_err = std::fmax(max_coef_err, max_coef_err_i);
	}

	return (max_rho_err <= tol_rho)&&(max_coef_err <= tol_coef);
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

	double tol_rho = 1e-16, tol_coef = 1e-04;
	test_by_tables4_9(tol_rho, tol_coef, true);

	std::system("pause");
	return 0;
}