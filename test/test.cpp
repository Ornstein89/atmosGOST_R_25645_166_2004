/* Test and example of use of the atmosGOST_R_25645_166_2004() function */

#include <iostream>
#include <string>
#include <chrono>

// определения для тестирования
#include <tables_4_9_GOST_R_25645_166_2004.h>
#include <tables_10_11_GOST_R_25645_166_2004.h>

// #define ATMOS_GOST_PARAMETERS_CHECK // включение проверки входных параметров функций, снижает производительность но облегчает поиск ошибок
// #define ATMOS_GOST_TEST // включение вывода тестовых данных из функции atmosGOST_R_25645_166_2004(), необходимо для тестирования, снижает производительность

#include <atmosGOST_R_25645_166_2004.h> // подключение

double F0_arr_test[7] = { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };

double rand_f(double a, double b) {
	return ((double)rand() / RAND_MAX) * (b - a) + a;
}

bool test_by_tables4_9(
	double const tol_rho,
	double const tol_coef,
	bool const verbose)
{
	if (verbose) {
		std::cout << std::endl << "=== Testing by matching GOST R 25645.166-2004 Tables 4-9 ==="
			<< std::endl;
		std::cout << "Result" << "\th,km" << "\t\tF107"
			<< "\t\trho erros"
			<< "\t\t\t\tK0-K4' erros"
			<< std::endl;
	}
	double max_rho_err = 0.0;
	double max_coef_err = 0.0;
	double max_rho_err_i = 0.0;
	double max_coef_err_i = 0.0;
	for (int i = 0; i < 7; i++) // перечисление по шкале F
	{
		for (int j = 0; j < 70; j++) // перечисление по шкале высоты
		{
			max_rho_err_i = 0.0;
			max_coef_err_i = 0.0;
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
			double err_rhoK1 = abs(test_glob_K1 - table_K1[j][i]);
			double err_rhoK2 = abs(test_glob_K2 - table_K2[j][i]);
			double err_rhoK3 = abs(test_glob_K3 - table_K3[j][i]);
			double err_rhoK4_1 = abs(test_glob_K4_1 - table_K4[j][i]);

			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK0);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK1);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK2);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK3);
			max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK4_1);

			if (verbose) {
				std::string rslt =
					((max_rho_err_i < tol_rho)&&(max_coef_err_i < tol_coef))
					? ("OK")
					: ("**ERR");
				std::cout << rslt << "\th = " << h
					<< ";\tF107 = " << F0_arr_test[i]
					/*<< ";\tmax_rho_error = " << max_rho_err_i
					<< ";\t\tmax_coef_error = " << max_coef_err_i
					<< std::endl*/
					<< "\t" /*<< table_rho_night[j][i] << "/"*/ << test_glob_rho_night
					<< "\t" /*<< table_K0[j][i] << "/"*/ << test_glob_K0
					<< "\t" /*<< table_K1[j][i] << "/"*/ << test_glob_K1
					<< "\t" /*<< table_K2[j][i] << "/"*/ << test_glob_K2
					<< "\t" /*<< table_K3[j][i] << "/"*/ << test_glob_K3
					<< "\t" /*<< table_K4[j][i] << "/"*/ << test_glob_K4_1
					<< std::endl;
			}
			max_rho_err = std::fmax(max_rho_err, max_rho_err_i);
			max_coef_err = std::fmax(max_coef_err, max_coef_err_i);
		}
	}

	return (max_rho_err <= tol_rho)&&(max_coef_err <= tol_coef);
}

bool test_by_tables10_11(
	double const tol_coef,
	bool const verbose)
{
	if (verbose) {
		std::cout << std::endl << "=== Testing by matching GOST R 25645.166-2004 Table 10 ==="
			<< std::endl;
		std::cout << "Result" << "\th,km" << "\t\tF107"
			<< "\t\tK4'' error" << std::endl;
	}
	double max_coef_err = 0.0, err_rhoK4_2 = 0.0;
	for (int i = 0; i < 7; i++) // перечисление по шкале F
	{
		for (int j = 0; j < 22; j++) // перечисление по шкале Kp
		{
			double Kp = Kp0 + j * Kp_step;
			double h = rand_f(120.0, 1500.0);
			double X0[3] = { 0.0 };
			X0[0] = h;

			double rho = atmosGOST_R_25645_166_2004(
				h,
				F0_arr_test[i],
				Kp,
				F0_arr_test[i], // 
				0.0,
				X0,
				0.0,
				0.0,
				0.0,
				0.0);

			err_rhoK4_2 = abs(test_glob_K4_2 - table_K4_2_24h[j][i]);

			if (verbose) {
				std::string rslt =
					(err_rhoK4_2 < tol_coef)
					? ("OK")
					: ("**ERR");
				std::cout << rslt
					<< "\th = " << h
					<< ";\tF107 = " << F0_arr_test[i]
					<< ";\terr_rhoK4_2 = " << err_rhoK4_2
					<< std::endl;
			}
		}
		max_coef_err = std::fmax(max_coef_err, err_rhoK4_2);
	}
	if (verbose) {
		std::cout << "\tMaximal K4'' error = " << max_coef_err << std::endl;
	}
	return (max_coef_err <= tol_coef);
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

	double tol_rho = 1e-10, tol_coef = 1.0e-03;
	test_by_tables4_9(tol_rho, tol_coef, true);
	test_by_tables10_11(tol_coef, true);

	/*
	for (int i = 0; i < 100; i++) {
		double Kp = 3.0;
		double h = rand_f(120.0, 1500.0);
		double F = rand_f(65.0, 260.0);
		double X0[3] = { 0.0 };
		X0[0] = h;
		std::cout << "Parameter F = " << F << "\t";
		double rho = atmosGOST_R_25645_166_2004(h, F, Kp, F, 0.0, X0, 0.0, 0.0, 0.0, 0.0);
		std::cout << std::endl;
	}*/

	std::system("pause");
	return 0;
}
