#pragma once

#ifndef ATMOS_GOST_H
#define ATMOS_GOST_H

double const H0 = 120.0;	//TODO	// высота начала шкалы таблиц 4-9 //TODO
double const H_end = 1500.0;	//TODO	// конечная высота шкалы таблиц 4-9 //TODO
double const H_step = 20.0;		// шаг шкалы таблиц 4-9 по высоте


/* Шкала опорных значений индекса солнечной активности F0,
 * по которой выстроены колонки таблиц 2-3 */
extern double F0_arr[7]; // { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };

double ApToKpLinearInterp(const double Ap);

double KpToApLinearInterp(const double Kp);

double calcF81(double const F107[81]);

int getClosestF0(double F81);

double rho_night(const double h_km, const double h_km_powers[7], const int n_col);

double K0_prime(const double h_km, const double h_km_powers[7], int n_col);

double K1_prime(const double h_km, const double h_km_powers[7], int n_col);

double K2_prime(const double h_km, const double h_km_powers[7], int n_col);

double K3_prime(const double h_km, const double h_km_powers[7], int n_col);

double K4_prime(const double h_km, const double h_km_powers[7], int n_col);

double K4_prime2_24h(const double Kp, int n_col);

double atmosGOST_R_25645_166_2004(
    double const h_km,
    double const F107,
    double const Kp,
    double const F81,
    double const DoY,
    double const X[3],
    double t_s,
    double S_rad,
    double alpha_rad,
    double delta_rad);

#endif // ATMOS_GOST_H
