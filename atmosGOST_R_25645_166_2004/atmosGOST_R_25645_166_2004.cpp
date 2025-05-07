#include "atmosGOST_R_25645_166_2004.h"
#include "table_1_GOST_R_25645_166_2004.h"
#include "table_2_GOST_R_25645_166_2004.h"
#include "table_3_GOST_R_25645_166_2004.h"

#define CONST_RHO_0		1.58868e-08	 // кг / м ^ 3, плотность ночной атмосферы на высоте 120км
#define CONST_OMEGA_Z	7.292115e-05 // рад/с, угловая скорость вращения Земли

// TODO учесть время запаздывания
#define CONST_dt_f107	1.7	 // запаздывание по F107
#define CONST_dt_Kp		0.6	 // запаздывание по Kp
#define CONST_dt_kpp	0.25 // сут, запаздывание по kpp

// #define _MATH_H
#define _USE_MATH_DEFINES
#include <cmath>


double F0_arr[7] = { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };

// транспонированные таблицы 2 и 3, позволяющие ссылаться на колонки
// без копирования элементов
// double a_high_transposed;
// double a_low_transposed[7][8] = {{0.0}};
// double b_high_transposed;
// double b_low_transposed[7][6] = {{0.0}};
// double c_high_transposed;
// double c_low_transposed[7][6] = {{0.0}};
// double d_high_transposed;
// double d_low_transposed[7][6] = {{0.0}};
// double e_high_transposed;
// double e_low_transposed[10][7] = {{0.0}};
// double l_high_transposed;
// double l_low_transposed[7][6] = {{0.0}};

void init_library()
{
    // транспонирование таблиц коэффициентов, позволяющее
    // использовать ссылку на колонку, например a[7] = a_high_table[n_col]
    // и не производить поэлементное копирование

    // TODO
}

const double tableAp[28] = {
    0,      2,      3,      4,      5,      6,      7,
    9,      12,     15,     18,     22,     27,     32,
    39,     48,     56,     67,     80,     94,     111,
    132,    154,	179,	207,	236,	300,	400};

/**
 * @brief Функция для перевода планетарного среднесуточного индекса
 * геомагнитной активности Ap [нТл] в квазилогарифмический среднесуточный
 * индекс Kp [баллы] по таблице А.1 ГОСТ Р 25645.166-2004
 * @param ap
 * @return
 */
double ApToKpLinearInterp(const double Ap)
{
    //TODO обработка ошибок
    int i;
    for(i=0; i < 27; i++){
        if(tableAp[i] <= Ap && Ap <= tableAp[i+1])
            break;
        if(i==26)
            return -1;
    }
    double result = (i*1.0/3.0)
                    + 1.0/3.0 * (Ap - tableAp[i]) / (tableAp[i+1] - tableAp[i]);
    return result;
}


/**
 * @brief Функция для перевода квазилогарифмического среднесуточного
 * индекса Kp [баллы] в планетарный среднесуточный индекс
 * геомагнитной активности Ap [нТл] в  по таблице А.1 ГОСТ Р 25645.166-2004
 * @param Kp
 * @return
 */
double KpToApLinearInterp(const double Kp)
{
    //TODO обработка ошибок
    double step = 1.0/3.0;
    int i = std::floor(Kp / step);
    double nearestLowKp = i * step;
    double result = tableAp[i]
                    + (tableAp[i+1] - tableAp[i]) * (Kp - nearestLowKp) / step;
    return result;
}


/** значения F107 начиная с 80 суток перед текущей датой до текущей даты */
double calcF81(double const F107[81])
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


/** выбор номера колонки F0 для таблиц 2 и 3 */
int getClosestF0(double F81)
{
    /* проверка F81 по допустимому диапазону */
    if(F81<0.0)
        return -1;

    if(F81 < (75.0+100)/2) // значение ближайшее к первому значению F0
        return 0;
    else if(F81 >= 225) // значение, ближайшее к последнему значению F0
        return 6;
    else // внутри таблицы - поиск ближайшего
        return std::round(F81/25) - 3;
}


/**
 * @brief rho_night - функция для расчёта плотности ночной атмосферы по
 * ГОСТ Р 25645.166-2004 п. 5.4
 * @param h_km
 * @param h_km_powers - массив степеней высоты [1, h, h^2, ..., h^6]
 * с целью оптимизации рассчитан один раз вне функции
 * @param n_col
 * @return
 */
double rho_night(const double h_km, const double h_km_powers[7], const int n_col)
{

    double a[7] = { 0.0 };

    // ГОСТ п.5.5: два высотных дипазона
    if (h_km > a_high_table[0][n_col]) { // если ">" - тест по таблице 4 проходит, если ">=" - не проходит, в ГОСТ нет комментариев, считать границу включительно или нет
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

    double polynom = 0.0;
    for (int i = 0; i < 7; i++) // сборка полинома
    {
        polynom += a[i] * h_km_powers[i];
    }

    double rho_night = CONST_RHO_0 * std::exp(polynom);
    return rho_night;
}


/**
 * @brief K0_prime - функция расчёта K'0 из ГОСТ Р 25645.166-2004 п. 5.8,
 * величины, характеризующей вековое изменение плотности ночной атмосферы
 * в 11-летнем цикле солнечной активности
 * @param h_km
 * @param h_km_powers - массив степеней высоты [1, h, h^2, ..., h^6]
 * с целью оптимизации рассчитан один раз вне функции
 * @param n_col
 * @return
 */
double K0_prime(const double h_km, const double h_km_powers[7], int n_col)
{
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

    double result = 0.0;

    for (int i = 0; i < 5; i++) // сборка полинома
    {
        result += l[i] * h_km_powers[i];
    };

    return result;
}


/**
 * @brief K1_prime - функция расчёта K'1 из ГОСТ Р 25645.166-2004 п. 5.8,
 * величины, характеризующей амплитуду суточного эффекта
 * @param h_km
 * @param h_km_powers - массив степеней высоты [1, h, h^2, ..., h^6]
 * с целью оптимизации рассчитан один раз вне функции
 * @param n_col
 * @return
 */
double K1_prime(const double h_km, const double h_km_powers[7], int n_col)
{
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

    double result = 0.0;
    for (int i = 0; i < 5; i++) // сборка полинома
    {
        result += c[i] * h_km_powers[i];
    };
    return result;
}


/**
 * @brief K2_prime - функция расчёта K'2 из ГОСТ Р 25645.166-2004 п. 5.8,
 * величины, характеризующей влияние полугодового эффекта
 * @param h_km
 * @param h_km_powers - массив степеней высоты [1, h, h^2, ..., h^6]
 * с целью оптимизации рассчитан один раз вне функции
 * @param n_col
 * @return
 */
double K2_prime(const double h_km, const double h_km_powers[7], int n_col)
{
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

    double result = 0.0;
    for (int i = 0; i < 5; i++) // сборка полинома
    {
        result += d[i] * h_km_powers[i];
    };
    return result;
}


/**
 * @brief K3_prime - функция расчёта K'3 из ГОСТ Р 25645.166-2004 п. 5.8,
 * величины, характеризующей влияние радиоизлучения Солнца
 * @param h_km
 * @param h_km_powers - массив степеней высоты [1, h, h^2, ..., h^6]
 * с целью оптимизации рассчитан один раз вне функции
 * @param n_col
 * @return
 */
double K3_prime(const double h_km, const double h_km_powers[7], int n_col)
{
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

    double result = 0.0;
    for (int i = 0; i < 5; i++) // сборка полинома
    {
        result += b[i] * h_km_powers[i];
    };
    return result;
}


/**
 * @brief K4_prime - функция расчёта K'4 из ГОСТ Р 25645.166-2004 п. 5.8,
 * величины, характеризующей среднесуточное влияние геомагнитной возмущённости
 * @param h_km
 * @param h_km_powers - массив степеней высоты [1, h, h^2, ..., h^6]
 * с целью оптимизации рассчитан один раз вне функции
 * вне функции
 * @param n_col
 * @return
 */
double K4_prime(const double h_km, const double h_km_powers[7], int n_col)
{
    // TODO после тестирования - избавиться от этого отдельного массива для быстродействия
    double e0to4[5] = { 0.0 };  // коэффициенты e0-e4
    if (h_km > e_high_table[0][n_col]) {
        for (int i = 0; i < 5; i++)
        {
            e0to4[i] = e_high_table[i + 1][n_col];
        }
    }
    else {
        for (int i = 0; i < 5; i++)
        {
            e0to4[i] = e_low_table[i + 1][n_col];
        }
    }

    double result = 0.0;
    for (int i = 0; i < 5; i++) // сборка полинома
    {
        result += e0to4[i] * h_km_powers[i];
    };
    return result;
}


/**
 * @brief K4_prime2 - функция расчёта K"4 из ГОСТ Р 25645.166-2004 п. 5.8,
 * ещё одной величины, характеризующей среднесуточное влияние
 * геомагнитной возмущённости (не путать с K"4 на трёхчасовом интервале,
 * рассчитанным через коэффициенты et)
 * @param Kp
 * @param n_col
 * @return
 */
double K4_prime2_24h(const double Kp, int n_col)
{
    double e5to8[4] = { 0.0 }; // коэффициенты e5-e8

    // коэффициенты e5-e8 одинаковы в обеих таблицах 2 и 3
    // TODO после тестирования - избавиться от этого отдельного массива для быстродействия
    for (int i = 0; i < 4; i++)
    {
        e5to8[i] = e_high_table[i + 1 + 5][n_col];
    }

    double Kp_powers[4] = { 0.0 }; // массив степеней Kp: [1, Kp, Kp^2, Kp^3]
    Kp_powers[0] = 1.0;
    for (int i = 1; i < 4; i++)
    {
        Kp_powers[i] = Kp_powers[i - 1] * Kp;
    }

    double result = 0.0; // второй множитель K4, K"4 ГОСТ Р 25645.166-2004 п. 5.8
    for (int i = 0; i < 4; i++) // сборка полинома
    {
        result += e5to8[i] * Kp_powers[i];
    }
    return result;
}


/**
 * @brief
 * [EN] Function atmosGOST_R_25645_166_2004() calculates uper Earth atmosphere
 * density by russian GOST R 25645.166-2004 model.
 * Returns
 * [RU] Функция atmosGOST_R_25645_166_2004() для расчёта плотности верхней
 * земной атмосферы по модели из российского стандарта ГОСТ Р 25645.166-2004.
 *
 * @param h_km - altitude above Earth ellipsoid, km, 120km<=h_km<=1500km /
 * высота над уровнем земного эллипсоида, км, от 120 до 1500км включительно
 * @param F107 - F10.7cm solar emission index, solar flux units (s.f.u.) /
 * индекс солнечной активности, единицы
 * @param Kp - quasi-logarithmic planetary 24-hour mean geomagnetic index, units /
 * квазилогарифмический планетарный среднесуточный индекс геомагнитной
 * активности, баллы
 * @param F81 - averaged weighted F10.7 for previous 80 days + current day /
 * усреднённый за 81 сутки (80 предыдущих + 1 текущие) и взвешенный индекс
 * солнечной активности
 * @param DoY - number of day from the beginning of the year / номер суток от
 * начала года
 * @param X[3] - x, y, z - geocentric greenwich coordinates, km / гринвичские
 * геоцентрические координаты точи пространства, км
 * @param t_s - universal time, s / всемирное время, с
 * @param S_rad - sidereal midnight time, rad / звёздное время в гринвичскую
 * полночь, рад
 * @param alpha_rad - right ascention of the Sun, rad / прямое восхождение
 * Солнца, рад
 * @param delta_rad - declination of the Sun, rad / склонение Солнца, рад
 * @return Density kg/m^3 if succeeded, -1 if failed. / В случае успешного
 * расчёта возвращает плотность в кг/куб.м, или -1 в случае ошибки.
 */
double atmosGOST_R_25645_166_2004(
    double const h_km,
    double const F107,
    double const Kp,
    double const F81,
    double const DoY,
    double const X[3],
    double const t_s,
    double const S_rad,
    double const alpha_rad,
    double const delta_rad)
{
    if((h_km < 120.0) || (h_km > 1500.0))
        return -1;
    if ((DoY < 0.0) || (DoY > 366.0))
        return -1.0;
    if ((alpha_rad < 0.0) || (alpha_rad > 2*M_PI))
        return -1.0;
    if ((delta_rad < -M_PI) || (delta_rad > M_PI))
        return -1.0;

    /* 2. Предварительный расчёт степеней от h_km, DoY, Kp для последующего
     * многократного использования:
     * h_km_powers = {1, h, h^2, h^3, ..., h^6} */
    double h_km_powers[7];
    h_km_powers[0] = 1.0;
	for (int i = 1; i < 7; i++)
	{
        h_km_powers[i] = h_km_powers[i - 1] * h_km;
	}

    // DOY = Day of the Year, doy_powers = {1, DoY, DoY^2, DoY^3, ..., DoY^8}
    double doy_powers[9] = { 0.0 };
    doy_powers[0] = 1.0;
	for (int i = 1; i < 9; i++)
	{
        doy_powers[i] = doy_powers[i - 1] * DoY;
	}

	// 3. Выбор опорного значения индекса солнечной активности F0
	// и номера колонки для таблицы 2 или 3
    int n_col = getClosestF0(F81);
	
    // 4. Выбор коэффициентов из таблиц, учитывая высоту (таблицы 2 либо 3 ГОСТа)
	// в том порядке, в котором они приведены в таблицах 2-3

	double n[3] = { 0.0 };
    if (h_km > c_high_table[0][n_col]) { // высота та же, что и для коэффициента c
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

    // коэффициент модели, равный углу запаздывания максимума плотности
    // по отношению к максимуму освещённости, рад
	double phi1 = 0.0;
    if (h_km > c_high_table[0][n_col]) // высота та же, что и для коэффициента c
		phi1 = phi1_high_table[n_col];
    else
		phi1 = phi1_low_table[n_col];


    /* для 3-х часового индекса k_pp геомагнитной возмущённости
     * TODO в отдельную функцию

	double et[4] = { 0.0 };
    if (h_km > e_high_table[0][n_col]) { // высота та же, что и для коэффициента e
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
    } */

	// 5. Расчёт коэффициентов K0-K4
	// 5.1 коэффициент K0, учитывающий изменение плотности атмосферы, связанное с отклонением

    double K0 = K0_prime(h_km, h_km_powers, n_col)
                * (F81 - F0_arr[n_col]) / F0_arr[n_col] + 1.0; // коэффициент в связи с отклонением F81 от F0

	// 5.2 коффициент K1, учитывающий суточный эффект в распределении плотности
    // (c' * h_km_powers(1:5) ) * (cos_phi ^ ( n * h_km_powers(1:3) )) / 2; // суточный коэффициент распределения плотности

	double cos_power = 0.0;
	for (int i = 0; i < 3; i++) // сборка полинома
        cos_power += n[i] * h_km_powers[i];

    // расстояние от центра гривничской СК
    double r = sqrt(X[0]* X[0] + X[1]*X[1] + X[2]*X[2]);

	// разность между долготой, для которой рассчитывают плотность атмосферы
	// и долготой с максимальным значением плотности в её суточном распределении, рад
	double beta_rad = alpha_rad - S_rad - CONST_OMEGA_Z * t_s + phi1;

	double cos_phi = 1 / r * (X[2] * std::sin(delta_rad)
		+ cos(delta_rad) * (X[0] * std::cos(beta_rad) + X[1] * sin(beta_rad)));
    double cos_half_phi = std::sqrt(1+cos_phi / 2);
    double K1 = K1_prime(h_km, h_km_powers, n_col)
                * std::pow(cos_half_phi, cos_power); // исправлена ошибка cos(phi)/2 -> cos(phi/2)

    // 5.3 Коэффициент K2, учитывающий полугодовой эффект
	double Ad = 0.0;
	for (int i = 0; i < 9; i++) // сборка полинома
        Ad += A_table[i] * doy_powers[i];

    double K2 = K2_prime(h_km, h_km_powers, n_col) * Ad;

    // 5.4 Коэффициент K3, учитывающий изменение плотности, связанное с отклонением F107 от F81
    double K3 = K3_prime(h_km, h_km_powers, n_col)
                * (F107 - F81) / (F81 + std::abs(F107 - F81));

    // 5.5 Коэффициент K4, учитывающий зависимость плотности атмосферы от
    // геомагнитной возмущенности при использовании среднесуточных (24h)
    // коэффициентов геомагнитной активности Kp
    double K4 = K4_prime(h_km, h_km_powers, n_col)
                * K4_prime2_24h(Kp, n_col);

    // Финальная формула
    return rho_night(h_km, h_km_powers, n_col) * K0 * (1 + K1 + K2 + K3 + K4);
}
