#pragma once

extern double test_glob_rho_night;
extern double test_glob_K0;
extern double test_glob_K1;
extern double test_glob_K2;
extern double test_glob_K3;
extern double test_glob_K4_1;
extern double test_glob_K4_2;

double atmosGOST_R_25645_166_2004(
	double const h_km,	// altitude above Earth ellipsoid 120<h_km<1500, km / высота над уровнем земного эллипсоида, от 120 до 1500, км 
	double const F107,	// F10.7 solar emission index / индекс солнечной активности
	double const Kp,	// квазилогарифмический планетарный среднесуточный индекс геомагнитной активности, баллы
	double const F81,	// averaged weighted F10.7 for previous 80 days + current day / усреднённый за 81 сутки (80 предыдущих + 1 текущие) и взвешенный индекс солнечной активности
	double const DoY,	// number of day from the beginning of the year / номер суток от начала года //TODO заменить на Mjd_TT
	double const X[3],	// x, y, z - geocentric greenwich coordinates, km / гринвичские координаты точи пространства, км
	double t_s,			// всемирное время, с
	double S_rad,		// sidereal midnight time, rad / звёздное время в гринвическую полночь, рад
	double alpha_rad,	// right ascention of the Sun, rad / прямое восхождение Солнца, рад
	double delta_rad);	// declination of the Sun, rad / склонение Солнца, рад

double calcF81(double const F107[81]);