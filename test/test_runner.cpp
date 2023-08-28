#include <gtest/gtest.h>

#include <vector>
#include <string>
// #define _MATH_H
// #define _USE_MATH_DEFINES
#include <cmath>

#include <atmosGOST_R_25645_166_2004.h>
#include <tables_4_9_GOST_R_25645_166_2004.h>
#include <tables_10_11_GOST_R_25645_166_2004.h>

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
            // относительная погрешность плотности - на 2 порядка меньше
            double tolerance = 0.5 * std::pow(10, std::floor(std::log10((long double)expected_density))-2);
            double density = rho_night(h_km, h_km_powers, F0i);
            EXPECT_TRUE(abs(density-expected_density) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected rho = " << expected_density
                << ", result = " << density;
        }
    }
}

/**
 *  Тест функции определения плотности ночной атмосферы K0_prime()
 *  (соответствует функции K0') по таблице 5 ГОСТ Р 25645.166-2004
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
            // относительная погрешность плотности - на 2 порядка меньше

            double tolerance = 0.001; //TODO

            EXPECT_TRUE(abs(value-expected_value) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << abs(expected_value-value);
        }
    }
}

/**
 *  Тест функции определения плотности ночной атмосферы K1_prime()
 *  (соответствует функции K1') по таблице 6 ГОСТ Р 25645.166-2004
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
            // относительная погрешность плотности - на 2 порядка меньше

            double tolerance = 0.001; //TODO


            EXPECT_TRUE(abs(value-expected_value) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции определения плотности ночной атмосферы K2_prime()
 *  (соответствует функции K2') по таблице 7 ГОСТ Р 25645.166-2004
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
            // относительная погрешность плотности - на 2 порядка меньше

            double tolerance = 0.001; //TODO

            EXPECT_TRUE(abs(value-expected_value) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции определения плотности ночной атмосферы K2_prime()
 *  (соответствует функции K3') по таблице 8 ГОСТ Р 25645.166-2004
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
            // относительная погрешность плотности - на 2 порядка меньше

            double tolerance = 0.001; //TODO

            EXPECT_TRUE(abs(value-expected_value) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции определения плотности ночной атмосферы K4_prime()
 *  (соответствует функции K4') по таблице 9 ГОСТ Р 25645.166-2004
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
            // относительная погрешность плотности - на 2 порядка меньше

            double tolerance = 0.001;   //TODO

            EXPECT_TRUE(abs(value-expected_value) < tolerance)
                << "Input: h_km = " << h_km << ", F0 = " << F0
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << abs(expected_value-value);
        }
    }
}


/**
 *  Тест функции определения плотности ночной атмосферы K4_prime2()
 *  (соответствует функции K4'') по таблице 10 ГОСТ Р 25645.166-2004
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
            // относительная погрешность плотности - на 2 порядка меньше

            double tolerance = 0.001; //TODO

            EXPECT_TRUE(abs(value-expected_value) < tolerance)
                << "Input: Kp = " << Kp
                << ", expected value = " << expected_value
                << ", value = " << value
                << ", error = " << abs(expected_value-value);
        }
    }
}

int main(int argc, char **argv) {

    // std::cout << "*** log10(1.62E-08) = " << floor(log10(1.62E-08));
    // std::cout << "*** log10(1.66E-13) = " << floor(log10(1.66E-13));
    // std::cout << "*** log10(0.0) = " << std::floor(std::log10(0.0));
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
