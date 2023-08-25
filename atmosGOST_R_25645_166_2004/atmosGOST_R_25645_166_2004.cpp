#include <atmosGOST_R_25645_166_2004.h>
#include <tables_1_3_GOST_R_25645_166_2004.h>

#define CONST_RHO_0		1.58868e-08		// кг / м ^ 3, плотность ночной атмосферы на высоте 120км
#define CONST_OMEGA_Z	7.292115e-05	// рад/с, угловая скорость вращения Земли
// TODO учесть время запаздывания
#define CONST_dt_f107	1.7		// запаздывание по F107
#define CONST_dt_Kp		0.6		// запаздывание по Kp
#define CONST_dt_kpp	0.25	// сут, запаздывание по kpp

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

double test_glob_rho_night;
double test_glob_K0;
double test_glob_K1;
double test_glob_K2;
double test_glob_K3;
double test_glob_K4_1;
double test_glob_K4_2;

double calcF81(double const F107[81]) //значения F107 начиная с 80 суток перед текущей датой до текущей даты
{
	/* формула по расчёту F81 из ГОСТ */
	double num = 0.0; // числитель
	double denom = 0.0; // знаменатель
	for (int i = -80; i < 1; i++)
	{
		double Wi = 1.0 + (0.5 * i) / 80.0;
		num += F107[i + 80] * Wi;
		denom += Wi;
	}
	double F81 = num / denom;
	return F81;
}

/*
	[EN] Function atmosGOST_R_25645_166_2004() calculates uper Earth atmosphere density
	by russian GOST R 25645.166-2004 model. Returns density kg/m^3 if succeeded, -1 if failed.
	[RU] Функция atmosGOST_R_25645_166_2004() для расчёта плотности верхней земной атмосферы
	по модели из российского стандарта ГОСТ Р 25645.166-2004. В случае успешного расчёта возвращает
	плотность в кг/куб.м, или -1 в случае ошибки.
*/
double atmosGOST_R_25645_166_2004(
	double const h_km,	// altitude above Earth ellipsoid 120<h_km<1500, km / высота над уровнем земного эллипсоида, от 120 до 1500, км 
	double const F107,	// F10.7 solar emission index / индекс солнечной активности
	double const Kp,	// квазилогарифмический планетарный среднесуточный индекс геомагнитной активности, баллы
	double const F81,	// averaged weighted F10.7 for previous 80 days + current day / усреднённый за 81 сутки (80 предыдущих + 1 текущие) и взвешенный индекс солнечной активности
	double const DoY,	// number of day from the beginning of the year / номер суток от начала года //TODO заменить на Mjd_TT
	double const X[3],	// x, y, z - geocentric greenwich coordinates, km / гринвичские координаты точи пространства, км
	double const t_s,	// всемирное время, с
	double const S_rad,	// sidereal midnight time, rad / звёздное время в гринвическую полночь, рад
	double const alpha_rad,	// right ascention of the Sun, rad / прямое восхождение Солнца, рад
	double const delta_rad)	// declination of the Sun, rad / склонение Солнца, рад
{
	#if defined(ATMOS_GOST_PARAMETERS_CHECK_ON)
	// 1. Проверка корректности ввода, незначительно замедляет производительность
	if ((h_km < 120.0) || (h_km > 1500.0)) {
		std::cerr << "Input parameters Altitude above ellipsoid h_km = " << h_km
			<<"km, while h must be 120 km <= h <= 1500 km." << std::endl;
		return -1.0;
	}
	if ((DoY < 0.0) || (DoY > 366.0)) {
		std::cerr << "Day of year DoY = " << DoY
			<< ", while h must be 0 <= DoY <= 366." << std::endl;
		return -1.0;
	}
	if ((alpha_rad < 0.0) || (alpha_rad > 2*M_PI)) {
		std::cerr << "Right ascention alpha_rad = " << alpha_rad
			<< ", while alpha_rad must be 0 <= alpha_rad <= 2*PI." << std::endl;
		return -1.0;
	}
	if ((delta_rad < -M_PI) || (delta_rad > M_PI)) {
		std::cerr << "Declination delta_rad = " << delta_rad
			<< ", while delta_rad must be -PI <= delta_rad <= PI." << std::endl;
		return -1.0;
	}
	#endif // PARAMETERS_CHECK_ON

	// 2. Предварительный расчёт степеней от h_km, DoY, Kp:
	// h_vec = {1, h, h^2, h^3, ..., h^6}
	double h_vec[7];
	h_vec[0] = 1.0;
	for (int i = 1; i < 7; i++)
	{
		h_vec[i] = h_vec[i - 1] * h_km;
	}

	// doy_vec = {1, DoY, DoY^2, DoY^3, ..., DoY^8}
	double doy_vec[9] = { 0.0 };
	doy_vec[0] = 1.0;
	for (int i = 1; i < 9; i++)
	{
		doy_vec[i] = doy_vec[i - 1] * DoY;
	}

	double Kp_vec[4] = { 0.0 };
	Kp_vec[0] = 1.0;
	for (int i = 1; i < 4; i++)
	{
		Kp_vec[i] = Kp_vec[i - 1] * Kp;
	}

	// 3. Выбор опорного значения индекса солнечной активности F0
	// и номера колонки для таблицы 2 или 3
	/*
	// альтернативный алгоритм поиска ближайшего  F0, показавший себя на 20-50% медленнее
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
	} */

	int n_col = 0;
	double diff1 = abs(F0_arr[n_col] - F81);
	double diff2 = diff1;
	while (n_col < 6)
	{
		n_col++;
		diff2 = abs(F0_arr[n_col] - F81);
		if (diff1 <= diff2) {
			--n_col;
			break;
		}
		diff1 = diff2;
	}
	/*#ifdef ATMOS_GOST_TEST_ON 
	std::cout << "F0 = " << F0_arr[n_col];
	#endif*/
	
	// 4. Выбор коэффициентов из таблиц, учитывая высоту (2 либо 3 таблицы ГОСТа)
	// в том порядке, в котором они приведены в таблицах 2-3
	double a[7] = { 0.0 };
	if (h_km > a_high_table[0][n_col]) {
		for (int i = 0; i < 7; i++)
		{
			a[i] = a_high_table[i + 1][n_col];
		}
	}
	else {
		for (int i = 0; i < 7; i++)
		{
			a[i] = a_low_table[i + 1][n_col];
		}
	}

	double b[5] = { 0.0 };
	if (h_km > b_high_table[0][n_col]) {
		for (int i = 0; i < 5; i++)
		{
			b[i] = b_high_table[i + 1][n_col];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			b[i] = b_low_table[i + 1][n_col];
		}
	}

	double c[5] = { 0.0 };
	if (h_km > c_high_table[0][n_col]) {
		for (int i = 0; i < 5; i++)
		{
			c[i] = c_high_table[i + 1][n_col];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			c[i] = c_low_table[i + 1][n_col];
		}
	}

	double n[3] = { 0.0 };
	if (h_km > c_high_table[0][n_col]) { // высота именно из таблицы коэффициента c
		for (int i = 0; i < 3; i++)
		{
			n[i] = n_high_table[i][n_col];
		}
	}
	else {
		for (int i = 0; i < 3; i++)
		{
			n[i] = n_low_table[i][n_col];
		}
	}

	// коэффициент модели, равный углу запаздывания максимума плотности по отношению к максимуму освещённости, рад
	double phi1 = 0.0;
	if (h_km > c_high_table[0][n_col]) { // высота именно из таблицы коэффициента c
		phi1 = phi1_high_table[n_col];
	}
	else {
		phi1 = phi1_low_table[n_col];
	}

	double d[5] = { 0.0 };
	if (h_km > d_high_table[0][n_col]) {
		for (int i = 0; i < 5; i++)
		{
			d[i] = d_high_table[i + 1][n_col];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			d[i] = d_low_table[i + 1][n_col];
		}
	}

	double e[9] = { 0.0 };
	if (h_km > e_high_table[0][n_col]) {
		for (int i = 0; i < 9; i++)
		{
			e[i] = e_high_table[i + 1][n_col];
		}
	}
	else {
		for (int i = 0; i < 9; i++)
		{
			e[i] = e_low_table[i + 1][n_col];
		}
	}

	double et[4] = { 0.0 };
	if (h_km > e_high_table[0][n_col]) { // высота именно из таблицы коэффициента e
		for (int i = 0; i < 4; i++)
		{
			et[i] = et_high_table[i][n_col];
		}
	}
	else {
		for (int i = 0; i < 4; i++)
		{
			et[i] = et_low_table[i][n_col];
		}
	}

	double l[5] = { 0.0 };
	if (h_km > l_high_table[0][n_col]) {
		for (int i = 0; i < 5; i++)
		{
			l[i] = l_high_table[i + 1][n_col];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			l[i] = l_low_table[i + 1][n_col];
		}
	}

	// 5. Расчёт коэффициентов K0-K4
	// 5.1 коэффициент K0, учитывающий изменение плотности атмосферы, связанное с отклонением
	double K0 = 0.0;
	for (int i = 0; i < 5; i++) // сборка полинома
	{
		K0 += l[i] * h_vec[i];
	};
	#ifdef ATMOS_GOST_TEST_ON 
	test_glob_K0 = K0;
	#endif
	K0 *= (F81 - F0_arr[n_col]) / F0_arr[n_col] + 1.0; // коэффициент в связи с отклонением F81 от F0

	// 5.2 коффициент K1, учитывающий суточный эффект в распределении плотности
	double K1 = 0.0; // (c' * h_vec(1:5) ) * (cos_phi ^ ( n * h_vec(1:3) )) / 2; // суточный коэффициент распределения плотности
	for (int i = 0; i < 5; i++) // сборка полинома
	{
		K1 += c[i] * h_vec[i];
	};
	#ifdef ATMOS_GOST_TEST_ON 
	test_glob_K1 = K1;
	#endif
	double cos_power = 0.0;
	for (int i = 0; i < 3; i++) // сборка полинома
	{
		cos_power += n[i] * h_vec[i];
	};

	double r = sqrt(X[0]* X[0] + X[1]*X[1] + X[2]*X[2]); // расстояние от центра гривничской СК

	// разность между долготой, для которой рассчитывают плотность атмосферы
	// и долготой с максимальным значением плотности в её суточном распределении, рад
	double beta_rad = alpha_rad - S_rad - CONST_OMEGA_Z * t_s + phi1;

	double cos_phi = 1 / r * (X[2] * sin(delta_rad)
		+ cos(delta_rad) * (X[0] * cos(beta_rad) + X[1] * sin(beta_rad)));
	K1 *= pow(cos_phi, cos_power)/2;

	// 5.3 коэффициент K2, учитывающий полугодовой эффект
	double K2 = 0.0;
	for (int i = 0; i < 5; i++) // сборка полинома
	{
		K2 += d[i] * h_vec[i];
	}
	#ifdef ATMOS_GOST_TEST_ON 
	test_glob_K2 = K2;
	#endif
	double Ad = 0.0;
	for (int i = 0; i < 9; i++) // сборка полинома
	{
		Ad += A[i] * doy_vec[i]; // полугодовой коэффициент
	}
	K2 *= Ad;

	// 5.4 коэффициент K3, учитывающий изменение плотности, связанное с отклонением F107 от F81
	double K3 = 0.0;
	for (int i = 0; i < 5; i++) // сборка полинома
	{
		K3 += b[i] * h_vec[i];
	}
	#ifdef ATMOS_GOST_TEST_ON 
	test_glob_K3 = K3;
	#endif

	K3 *= (F107 - F81) / (F81 + abs(F107 - F81)); // изменение плотности в связи с отклонением F10.7 от F81

	// 5.5 Расчёт коэффициента K4, учитывающий зависимость плотности атмосферы от геомагнитной возмущенности
	// при использовании среднесуточных коэффициентов геомагнитной активности
	// e[5]...e[8]  и  Kp, при использовании 3-х часовых - et[5]...et[8] и kpp
	double K4_1 = 0.0; // первый множитель K4
	for (int i = 0; i < 5; i++) // сборка полинома
	{
		K4_1 += e[i] * h_vec[i];
	}

	double K4_2 = 0.0; // второй множитель K4
	for (int i = 0; i < 4; i++) // сборка полинома
	{
		K4_2 += e[i+5] * Kp_vec[i];
	}
	#ifdef ATMOS_GOST_TEST_ON 
	test_glob_K4_1 = K4_1;
	test_glob_K4_2 = K4_2;
	#endif
	double K4 = K4_1 * K4_2;

	// 5.6 финальная формула
	double polynom = 0.0;
	for (int i = 0; i < 7; i++) // сборка полинома
	{
		polynom += a[i] * h_vec[i];
	}
	double rho_night = CONST_RHO_0 * exp(polynom);

	#ifdef ATMOS_GOST_TEST_ON 
	test_glob_rho_night = rho_night;
	#endif

	return rho_night * K0 * (1 + K1 + K2 + K3 + K4);
}
