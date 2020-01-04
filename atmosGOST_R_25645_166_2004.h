#pragma once
double atmosGOST_R_25645_166_2004(double const h_km, // altitude, km; 120<h_km<1500
	double const F107,
	double const F81, // 
	double const DoY,
	double const X[3], // x, y, z - гринвичские координаты точи пространства, км
	double t, // всемирное время, с
	double S, // звёздное время в гринвическую полночь, рад
	double alpha, // прямое восхождение Солнца, рад
	double delta // склонение Солнца, рад
);

double calcF81(double const F107[81]);