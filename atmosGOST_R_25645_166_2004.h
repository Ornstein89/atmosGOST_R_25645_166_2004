#pragma once

extern double test_glob_rho_night;
extern double test_glob_K0;
extern double test_glob_K1;
extern double test_glob_K2;
extern double test_glob_K3;
extern double test_glob_K4_1;
extern double test_glob_K4_2;

double atmosGOST_R_25645_166_2004(
	double const h_km,	// altitude above Earth ellipsoid 120<h_km<1500, km / ������ ��� ������� ������� ����������, �� 120 �� 1500, �� 
	double const F107,	// F10.7 solar emission index / ������ ��������� ����������
	double const Kp,	// �������������������� ����������� �������������� ������ ������������ ����������, �����
	double const F81,	// averaged weighted F10.7 for previous 80 days + current day / ���������� �� 81 ����� (80 ���������� + 1 �������) � ���������� ������ ��������� ����������
	double const DoY,	// number of day from the beginning of the year / ����� ����� �� ������ ���� //TODO �������� �� Mjd_TT
	double const X[3],	// x, y, z - geocentric greenwich coordinates, km / ����������� ���������� ���� ������������, ��
	double t_s,			// ��������� �����, �
	double S_rad,		// sidereal midnight time, rad / ������� ����� � ������������ �������, ���
	double alpha_rad,	// right ascention of the Sun, rad / ������ ����������� ������, ���
	double delta_rad);	// declination of the Sun, rad / ��������� ������, ���

double calcF81(double const F107[81]);