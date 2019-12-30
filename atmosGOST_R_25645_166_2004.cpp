#include "atmosGOST_R_25645_166_2004.h"
#include <iostream>

// �������
double const a_high_table[8][7] = {
	{ 500,			500,			500,			500,			500,			500,			500 }, // ������ ������� ��������� ���������, ��
	{ 17.8781,		-2.54909,		-13.9599,		-23.3079,		-14.7264,		-4.912,			-5.40952 }, // a1
	{ -0.132025,	0.0140064,		0.0844951,		0.135141,		0.0713256,		0.0108326,		0.00550749 }, // a2
	{ 2.27717E-04,	-1.6946E-04,	-3.28875E-04,	-4.20802E-04,	-2.28015E-04,	-8.10546E-05,	-3.78851E-05 }, // a3
	{ -2.2543E-07,	3.27196E-07,	5.05918E-07,	5.73717E-07,	2.8487E-07,		1.15712E-07,	2.4808E-08 }, // a4
	{ 1.33574E-10,	-2.8763E-10,	-3.92299E-10,	-4.03238E-10,	-1.74383E-10,	-8.13296E-11,	4.92183E-12 }, // a5
	{ -4.50458E-14,	1.22625E-13,	1.52279E-13,	1.42846E-13,	5.08071E-14,	3.04913E-14,	-8.65011E-15 }, // a6
	{ 6.72086E-18,	-2.05736E-17,	-2.35576E-17,	-2.01726E-17,	-5.34955E-18,	-4.94989E-18,	1.9849E-18 } }; // a7

double const a_low_table[8][7] = {
	{ 120,	120,	120,	120,	120,	120,	120 },
	{ 26.8629,	27.4598,	28.6395,	29.6418,	30.1671,	29.7578,	30.7854 },
	{ -0.451674,	-0.463668,	-0.490987,	-0.514957,	-0.527837,	-0.517915,	-0.545695 },
	{ 2.90397E-03,	2.974E-03,	3.20649E-03,	3.41926E-03,	3.53211E-03,	3.42699E-03,	3.70328E-03 },
	{ -1.06953E-05,	-1.0753E-05,	-1.1681E-05,	-1.25785E-05,	-1.30227E-05,	-1.24137E-05,	-1.37072E-05 },
	{ 2.21598E-08,	2.17059E-08,	2.36847E-08,	2.5727E-08,	2.66455E-08,	2.48209E-08,	2.80614E-08 },
	{ -2.42941E-11,	-2.30249E-11,	-2.51809E-11,	-2.75874E-11,	-2.85432E-11,	-2.58413E-11,	-3.00184E-11 },
	{ 1.09926E-14,	1.00123E-14,	1.09536E-14,	1.21091E-14,	1.25009E-14,	1.09383E-14,	1.31142E-14 } };

double const b_high_table[6][7] = {
	{ 600,	660,	760,	800,	860,	900,	1000 },
	{ 23.1584,	33.2732,	39.1961,	43.2469,	49.5738,	11.278,	-52.6184 },
	{ -0.0802147,	-0.111099,	-0.12352,	-0.126973,	-0.138613,	0.00143478,	0.214689 },
	{ 1.05824E-04,	1.41421E-04,	1.49015E-04,	1.42637E-04,	1.47851E-04,	-3.69846E-05,	-2.94882E-04 },
	{ -6.15036E-08,	-7.94952E-08,	-7.9705E-08,	-7.09985E-08,	-6.96361E-08,	3.58318E-08,	1.71171E-07 },
	{ 1.32453E-11,	1.65836E-11,	1.58772E-11,	1.31646E-11,	1.21595E-11,	-9.91225E-12,	-3.60582E-11 } };

double const b_low_table[6][7] = {
	{ 120,	120,	120,	120,	120,	120,	120 },
	{ 0.0687894,	0.15073,	0.0479451,	0.0223448,	-0.00326391,	-0.0514749,	-0.107255 },
	{ -0.00284077,	-0.00400889,	-0.00239453,	-0.0019798,	-0.00159869,	-0.000921059,	-0.000174343 },
	{ 1.83922E-05,	2.43937E-05,	1.70335E-05,	1.54101E-05,	1.40443E-05,	1.15147E-05,	9.02759E-06 },
	{ 9.19605E-09,	-9.92772E-09,	-1.31626E-09,	-2.3543E-09,	-3.02287E-09,	-1.22901E-09,	-3.16512E-10 },
	{ -4.16873E-11,	-1.82239E-11,	-1.74032E-11,	-1.24994E-11,	-9.2016E-12,	-8.13104E-12,	-6.14E-12 } };

double const c_high_table[6][7] = {
	{ 640,	700,	760,	820,	860,	920,	980 },
	{ 50.5034,	61.624,	53.2623,	18.2236,	-31.8442,	-48.7208,	-147.859 },
	{ -0.170541,	-0.192967,	-0.144342,	-0.00840024,	0.168327,	0.222996,	0.531652 },
	{ 2.17232E-04,	2.28061E-04,	1.4659E-04,	-3.88E-05,	-2.62603E-04,	-3.21884E-04,	-6.71937E-04 },
	{ -1.21902E-07,	-1.18715E-07,	-6.46443E-08,	4.31384E-08,	1.65454E-07,	1.91495E-07,	3.64787E-07 },
	{ 2.54037E-11,	2.29638E-11,	1.04227E-11,	-1.23832E-11,	-3.69355E-11,	-4.08067E-11,	-7.26268E-11 } };

double const c_low_table[6][7] = {
	{ 120,	120,	120,	120,	120,	120,	120 },
	{ -1.04825,	-0.93106,	-0.820867,	-0.744047,	-0.722471,	-0.687482,	-0.739984 },
	{ 0.0166305,	0.0141537,	0.0119916,	0.0104743,	0.00980317,	0.00916594,	0.00952854 },
	{ -9.24263E-05,	-7.29862E-05,	-5.79835E-05,	-4.78544E-05,	-4.25245E-05,	-3.80932E-05,	-3.62727E-05 },
	{ 2.72382E-07,	2.00294E-07,	1.50707E-07,	1.18513E-07,	9.95544E-08,	8.51275E-08,	7.3887E-08 },
	{ -2.41355E-10,	-1.62006E-10,	-1.13026E-10,	-8.31498E-11,	-6.55175E-11,	-5.29972E-11,	-4.23907E-11 } };

double const d_high_table[6][7] = {
	{ 1500,	1500,	1500,	1500,	1500,	1500,	1500 },
	{ -0.351899,	-0.047813,	0.20981,	0.265174,	0.23047,	0.170074,	0.088141 },
	{ 0.00577056,	0.00380813,	0.00262881,	0.00275836,	0.00338331,	0.00406131,	0.00468253 },
	{ 9.95819E-07,	4.22771E-06,	4.24379E-06,	2.08668E-06,	-5.52305E-07,	-2.82114E-06,	-4.24609E-06 },
	{ -7.25324E-09,	-8.66826E-09,	-6.67328E-09,	-3.69543E-09,	-8.23607E-10,	1.38369E-09,	2.53509E-09 },
	{ 2.9759E-12,	3.06712E-12,	2.13496E-12,	1.11862E-12,	2.21349E-13,	-4.27908E-13,	-7.29031E-13 } };

double const d_low_table[6][7] = {
	{ 120,	120,	120,	120,	120,	120,	120 },
	{ -0.351899,	-0.047813,	0.20981,	0.265174,	0.23047,	0.170074,	0.088141 },
	{ 0.00577056,	0.00380813,	0.00262881,	0.00275836,	0.00338331,	0.00406131,	0.00468253 },
	{ 9.95819E-07,	4.22771E-06,	4.24379E-06,	2.08668E-06,	-5.52305E-07,	-2.82114E-06,	-4.24609E-06 },
	{ -7.25324E-09,	-8.66826E-09,	-6.67328E-09,	-3.69543E-09,	-8.23607E-10,	1.38369E-09,	2.53509E-09 },
	{ 2.9759E-12,	3.06712E-12,	2.13496E-12,	1.11862E-12,	2.21349E-13,	-4.27908E-13,	-7.29031E-13 } };

double const l_high_table[6][7] = {
	{ 640,	660,	740,	800,	860,	900,	900 },
	{ 48.6536,	54.4867,	60.1267,	47.0996,	50.6174,	8.01942,	-15.5728 },
	{ -1.70291E-01,	-1.78298E-01,	-1.83144E-01,	-1.2526E-01,	-1.29047E-01,	1.85302E-02,	9.36704E-02 },
	{ 2.26242E-04,	2.22725E-04,	2.12481E-04,	1.26352E-04,	1.24842E-04,	-6.14733E-05,	-1.49036E-04 },
	{ -1.32032E-07,	-1.227E-07,	-1.08497E-07,	-5.51584E-08,	-5.24993E-08,	4.97674E-08,	9.42151E-08 },
	{ 2.85193E-11,	2.51316E-11,	2.0571E-11,	8.75272E-12,	8.08272E-12,	-1.26162E-11,	-2.0961E-11 } };

double const l_low_table[6][7] = {
	{ 120,	120,	120,	120,	120,	120,	120 },
	{ -0.407768,	-0.902739,	-0.733037,	-1.31444,	-1.20026,	-1.52158,	-1.67664 },
	{ 0.00148506,	0.00826803,	0.00523396,	0.0133124,	0.0114087,	0.015704,	0.0177194 },
	{ 1.25357E-05,	-1.25448E-05,	6.35667E-06,	-2.55585E-05,	-1.47324E-05,	-3.02859E-05,	-3.69498E-05 },
	{ 3.77311E-08,	6.12853E-08,	1.09065E-08,	5.43981E-08,	2.7804E-08,	4.57668E-08,	5.09134E-08 },
	{ -7.78953E-11,	-7.07966E-11,	-2.61427E-11,	-4.33784E-11,	-2.2632E-11,	-2.82926E-11,	-2.82878E-11 } };

double n_high_table[4][7] = {
	{ 640,	700,	760,	820,	860,	920,	980 },
	{ 2.058,	2.058,	2.058,	2.058,	2.058,	2.058,	2.058 },
	{ 5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03 },
	{ -4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06 } };

double const n_low_table[4][7] = {
	{ 120,	120,	120,	120,	120,	120,	120 },
	{ 2.058,	2.058,	2.058,	2.058,	2.058,	2.058,	2.058 },
	{ 5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03,	5.887E-03 },
	{ -4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06,	-4.012E-06 } };

double const A[9] = { -2.53418e-02, -2.44075e-03, 3.08389e-06,
2.90115e-06, -4.99606e-08, 3.36327e-10,
-1.0966e-12, 1.73227e-15, -1.06271e-18 };

double atmosGOST_R_25645_166_2004(double const h, // � - ������
	double const F107,
	double const F81,
	double const doy) // ���� �� ������ ����, �������� ��� ��������� ���������� �������, ����� ������ ���� ������� �� Mjd_TT
// �������� �������� - rho, ��/�^3 - ���������
{
	// �������������
	double hkm = h / 1000;

	if (hkm < 120)
	{
		std::cerr << "h < 120 km" << std::endl;
		return NAN;
	}
	else if (hkm > 1500)
	{
		std::cerr << "h > 1500 km" << std::endl;
		return NAN;
	}

	double h_vec[7];
	h_vec[0] = 1.0;
	for (int i = 1; i < 7; i++)
	{
		h_vec[i] = h_vec[i - 1] * hkm;
	}

	double doy_vec[9] = { 0.0 };
	doy_vec[0] = 1.0;
	for (int i = 1; i < 9; i++)
	{
		doy_vec[i] = doy_vec[i - 1] * doy;
	}

	double F0 = F81;
	if (F81 > 250) F0 = 250;
	else if (F81 < 75) F0 = 75;
	else F0 = round(F81 / 25) * 25;

	double rho0 = 1.58868e-08; // �� / � ^ 3, ��������� ������ ��������� �� ������ 120��
	/*
	//TODO
	{
		Kp = ;
	Kp_vec = [1 Kp Kp ^ 2 Kp ^ 3];

	alpha = ;
	S = ;
	omega = ;
	t = ;
	phi1 = ;

	beta = alpha - S - omega * t + phi1;

	cos_phi = 1 / r * (z * sin(delta) + cos(delta) * (x * cos(beta) + y * sin(beta)));
	*/
	// ����� ������������� �� ������
	int F_index = (F0 - 50) / 25;

	double a[7] = { 0.0 };
	if (hkm < a_high_table[0][F_index]) {
		for (int i = 0; i < 7; i++)
		{
			a[i] = a_low_table[i + 1][F_index];
		}
	}
	else {
		for (int i = 0; i < 7; i++)
		{
			a[i] = a_high_table[i + 1][F_index];
		}
	}

	double b[5] = { 0.0 };
	if (hkm < b_high_table[0][F_index]) {
		for (int i = 0; i < 5; i++)
		{
			b[i] = b_low_table[i + 1][F_index];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			b[i] = b_high_table[i + 1][F_index];
		}
	}

	double c[5] = { 0.0 };
	if (hkm < c_high_table[0][F_index]) {
		for (int i = 0; i < 5; i++)
		{
			c[i] = c_low_table[i + 1][F_index];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			c[i] = c_high_table[i + 1][F_index];
		}
	}

	double d[5] = { 0.0 };
	if (hkm < d_high_table[0][F_index]) {
		for (int i = 0; i < 5; i++)
		{
			d[i] = d_low_table[i + 1][F_index];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			d[i] = d_high_table[i + 1][F_index];
		}
	}

	double l[5] = { 0.0 };
	if (hkm < l_high_table[0][F_index]) {
		for (int i = 0; i < 5; i++)
		{
			l[i] = l_low_table[i + 1][F_index];
		}
	}
	else {
		for (int i = 0; i < 5; i++)
		{
			l[i] = l_high_table[i + 1][F_index];
		}
	}

	double n[3] = { 0.0 };
	if (hkm < n_high_table[0][F_index]) {
		for (int i = 0; i < 3; i++)
		{
			n[i] = n_low_table[i + 1][F_index];
		}
	}
	else {
		for (int i = 0; i < 3; i++)
		{
			n[i] = n_high_table[i + 1][F_index];
		}
	}

	// ������ ���������
	// ���������� ������������� ��� �������� :
	// �������� ������� a = fliplr(a)
	// K0 = 1 + polyval(a, hkm);

	double K0 = 0.0;
	for (int i = 0; i < 5; i++)
	{
		K0 += l[i] * h_vec[i];
	};
	K0 = K0 * (F81 - F0) / F0 + 1.0; // ����������� � ����� � ����������� F81 �� F0

	double K1 = 0.0; // (c' * h_vec(1:5) ) * (cos_phi ^ ( n * h_vec(1:3) )) / 2; // �������� ����������� ������������� ���������

	double K2 = 0.0;
	for (int i = 0; i < 5; i++)
	{
		K2 += d[i] * h_vec[i];
	}
	double Ad = 0.0;
	for (int i = 0; i < 9; i++)
	{
		Ad += A[i] * doy_vec[i]; // ����������� �����������
	}
	K2 *= Ad;

	double K3 = 0.0;
	for (int i = 0; i < 5; i++)
	{
		K3 += b[i] * h_vec[i];
	}

	K3 *= (F107 - F81) / (F81 + abs(F107 - F81)); // ��������� ��������� � ����� � ����������� F10.7 �� F81

	//TODO
	double K4 = 0.0; // (e(1:5) * h_vec(1:5)) * (e(6:9) * Kp_vec); % ����������� �� ������������ �������������

	double pwr = 0.0;
	for (int i = 0; i < 7; i++)
	{
		pwr += a[i] * h_vec[i];
	}
	double rho_night = rho0 * exp(pwr);

	double rho = rho_night * K0 * (1 + K1 + K2 + K3 + K4);

	return rho;
}