#include <gtest/gtest.h>

#include <vector>
#include <string>
#include <cmath>

#include <atmosGOST_R_25645_166_2004.h>
#include <tables_4_9_GOST_R_25645_166_2004.h>
#include <tables_10_11_GOST_R_25645_166_2004.h>

// таблица построена по scipy.interpolate.interp1d
const double testTableApKp[2][30] = {
    /*Ap:*/ {	10.37, 36.06, 38.94, 40.20, 44.23, 58.17, 59.05, 102.45, 110.48, 111.97, 118.75, 141.43, 142.72, 157.54, 167.12, 175.29, 176.17, 183.33, 185.29, 197.83, 200.95, 209.42, 232.55, 243.41, 276.94, 299.42, 331.93, 363.72, 385.63, 389.33	},
    /*Kp:*/ {	2.4856, 4.5267, 4.6638, 4.7111, 4.8604, 5.3991, 5.4258, 6.4990, 6.6565, 6.6821, 6.7897, 7.1429, 7.1624, 7.3805, 7.5083, 7.6172, 7.6289, 7.7182, 7.7415, 7.8908, 7.9280, 8.0278, 8.2937, 8.3719, 8.5466, 8.6636, 8.7731, 8.8791, 8.9521, 8.9644	}
};

TEST(ApToKpLinearInterp, positive)
{
    for(int i = 0; i < 30; i++){
        EXPECT_NEAR(testTableApKp[1][i],
                    ApToKpLinearInterp(testTableApKp[0][i]),
                    0.00005);
    }
}

TEST(KpToApLinearInterp, positive)
{
    for(int i = 0; i < 30; i++){
        EXPECT_NEAR(testTableApKp[0][i],
                    KpToApLinearInterp(testTableApKp[1][i]),
                    0.005);
    }
}

/**
 * Тест функции определения номера колонки для таблиц 2 и 3 по значению F81
 */
TEST(getClosestF0, positive)
{
    struct getClosestF0_inout {
        double F81;
        int reuslt;
    };

    std::vector<getClosestF0_inout> in_out = {
        {-0.1, -1},
        {0.0, 0},{60.0, 0},{80.0, 0},
        {90.0, 1}, {110.0, 1},
        {115.0, 2}, {135.0, 2},
        {145.0, 3}, {155.0, 3},
        {170.0, 4}, {180.0, 4},
        {190.0, 5}, {210.0, 5},
        {225.0, 6}, {240.0, 6}, {260.0, 6},
    };

    // for(const getClosestF0_inout & testcase : in_out){
    for(int i = 0; i < in_out.size(); i++){
        getClosestF0_inout testcase = in_out[i];
        int tmp_result = getClosestF0(testcase.F81);
        EXPECT_TRUE(testcase.reuslt == tmp_result)
            << "Input: " << testcase.F81
            << ", expected: " << testcase.reuslt
            << ", result = " << tmp_result;
    }
}


/**
 *  Тест функции определения плотности ночной атмосферы rho_night()
 *  по таблице 4 ГОСТ Р 25645.166-2004
*/
TEST(rho_night_table4, positive)
{
    for(int hi = 0; hi < 70; hi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double h_km = H0 + hi * H_step;
            double F0 = F0_arr[F0i];


            double h_km_powers[7];
            h_km_powers[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_km_powers[i] = h_km_powers[i - 1] * h_km;
            }
            double expected_density = table_rho_night[hi][F0i];

            /* относительная погрешность плотности - погрешность округления
             * (0,5 от последнего значащего разряда) плотности из таблицы 4 */
            double tolerance = 0.5 * std::pow(10, std::floor(std::log10((long double)expected_density))-2);
            double density = rho_night(h_km, h_km_powers, F0i);
            EXPECT_TRUE(std::abs(density-expected_density) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected rho = " << expected_density
                << ", result = " << density;
        }
    }
}


/**
 * @brief compare_as_double - процедура сравнения рассчитанных значений коэффициентов
 * Kn' (K0'...K4'') с ожидаемыми значениями из таблиц 5-10, основанная на
 * сравнении разности чисел в формате double с ожидаемой погрешностью
 * округления коэффициентов в таблицах ГОСТ, равной 0,0005
 * @param ecpected - ожидаемое значение коэффициента Kn' (K0'...K4'') из таблиц
 * @param calculated - рассчитанное значение
 * @return true если величины различаются менее, чем на величину погрешности
 */
bool compare_as_double(double expected, double calculated)
{
    /* относительная погрешность - погрешность округления
     * 0,0005 (т.е. 0,5 от последнего значащего разряда)
     * Kn' из таблиц 5-10 - ТЕСТ НЕ ПРОХОДИТ кроме K3';
     * при 0,001 - не проходит только K2' */
    double tolerance = 0.0005;
    return std::abs(calculated-expected) < tolerance;
}


/**
 * @brief compare_as_fixed - процедура сравнения рассчитанных значений коэффициентов
 * Kn' (K0'...K4'') с ожидаемыми значениями из таблиц 4-10, основанная на
 * сравнении чисел, приведённых к фиксированной точности до 3 разряда
 * @param ecpected - ожидаемое значение коэффициента Kn' (K0'...K4'') из таблиц
 * @param calculated - рассчитанное значение
 * @return true если величины различаются менее, чем на величину погрешности
 */
bool compare_as_fixed(double expected, double calculated)
{
    int expected_int = static_cast<int>(std::round(1000.0 * expected));
    int calculated_int = static_cast<int>(std::round(1000.0 * calculated));
    return (expected_int==calculated_int);
}


/**
 *  Тест функции K0_prime() (или K0') векового изменения солнечной активности
 *  в 11-летнем цикле, по таблице 5 ГОСТ Р 25645.166-2004
*/
TEST(K0_prime_table5, positive)
{
    for(int hi = 0; hi < 70; hi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double h_km = H0 + hi * H_step;
            double F0 = F0_arr[F0i];

            double h_km_powers[7];
            h_km_powers[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_km_powers[i] = h_km_powers[i - 1] * h_km;
            }

            double expected_value = table_K0[hi][F0i];
            double value = K0_prime(h_km, h_km_powers, F0i);

            //EXPECT_TRUE(compare_as_double(expected_value, value))
            EXPECT_TRUE(compare_as_fixed(expected_value, value))
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << std::abs(expected_value-value);
        }
    }
}

/**
 *  Тест функции K1_prime() (K1') амплитуды солнечного эффекта
 *  по таблице 6 ГОСТ Р 25645.166-2004
*/
TEST(K1_prime_table6, positive)
{
    for(int hi = 0; hi < 70; hi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double h_km = H0 + hi * H_step;
            double F0 = F0_arr[F0i];

            double h_km_powers[7];
            h_km_powers[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_km_powers[i] = h_km_powers[i - 1] * h_km;
            }

            double expected_value = table_K1[hi][F0i];
            double value = K1_prime(h_km, h_km_powers, F0i);

            //EXPECT_TRUE(compare_as_double(expected_value, value))
            EXPECT_TRUE(compare_as_fixed(expected_value, value))
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << std::abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции K2_prime() (K2') влияния полугодового эффекта
 *  по таблице 7 ГОСТ Р 25645.166-2004
*/
TEST(K2_prime_table7, positive)
{
    for(int hi = 0; hi < 70; hi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double h_km = H0 + hi * H_step;
            double F0 = F0_arr[F0i];

            double h_km_powers[7];
            h_km_powers[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_km_powers[i] = h_km_powers[i - 1] * h_km;
            }

            double expected_value = table_K2[hi][F0i];
            double value = K2_prime(h_km, h_km_powers, F0i);

            //EXPECT_TRUE(compare_as_double(expected_value, value))
            EXPECT_TRUE(compare_as_fixed(expected_value, value))
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << std::abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции K3_prime() (K3') влияния радиоизлучения Солнца
 *  по таблице 8 ГОСТ Р 25645.166-2004
*/
TEST(K3_prime_table8, positive)
{
    for(int hi = 0; hi < 70; hi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double h_km = H0 + hi * H_step;
            double F0 = F0_arr[F0i];

            double h_km_powers[7];
            h_km_powers[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_km_powers[i] = h_km_powers[i - 1] * h_km;
            }

            double expected_value = table_K3[hi][F0i];
            double value = K3_prime(h_km, h_km_powers, F0i);

            //EXPECT_TRUE(compare_as_double(expected_value, value))
            EXPECT_TRUE(compare_as_fixed(expected_value, value))
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << std::abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции K4_prime() (K4') - первого множителя для расчёта влияния
 *  геомагнитной возмущённости, по таблице 9 ГОСТ Р 25645.166-2004
*/
TEST(K4_prime_table9, positive)
{
    for(int hi = 0; hi < 70; hi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double h_km = H0 + hi * H_step;
            double F0 = F0_arr[F0i];

            double h_km_powers[7];
            h_km_powers[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_km_powers[i] = h_km_powers[i - 1] * h_km;
            }

            double expected_value = table_K4[hi][F0i];
            double value = K4_prime(h_km, h_km_powers, F0i);

            //EXPECT_TRUE(compare_as_double(expected_value, value))
            EXPECT_TRUE(compare_as_fixed(expected_value, value))
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << std::abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции K4_prime2() (K4'') - второго множителя для расчёта влияния
 *  геомагнитной возмущённости, по таблице 10 ГОСТ Р 25645.166-2004
*/
TEST(K4_prime2_table10, positive)
{
    for(int Kpi = 0; Kpi < 22; Kpi++) // высоты
    {
        for(int F0i = 0; F0i < 7; F0i++)
        {
            double Kp = Kp0 + Kpi * Kp_step;
            double expected_value = table_K4_2_24h[Kpi][F0i];
            double value = K4_prime2_24h(Kp, F0i);

            //EXPECT_TRUE(compare_as_double(expected_value, value))
            EXPECT_TRUE(compare_as_fixed(expected_value, value))
                << "Input: Kp = " << Kp
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << std::abs(expected_value-value);
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
