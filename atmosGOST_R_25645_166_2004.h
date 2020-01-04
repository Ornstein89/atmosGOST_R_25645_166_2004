#pragma once
double atmosGOST_R_25645_166_2004(double const h_km, // altitude, km; 120<h_km<1500
	double const F107,
	double const F81, // 
	double const DoY,
	double const X[3], // x, y, z - ����������� ���������� ���� ������������, ��
	double t, // ��������� �����, �
	double S, // ������� ����� � ������������ �������, ���
	double alpha, // ������ ����������� ������, ���
	double delta // ��������� ������, ���
);

double calcF81(double const F107[81]);