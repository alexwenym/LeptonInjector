
#include <math.h>
#include <cmath>
#include <iostream>

#include <gtest/gtest.h>

#include "LeptonInjector/math/Polynomial.h"

using namespace LI::math;

double f1(double x){
    return pow(x, 3) - 4;
}

double df1(double x){
    return 3. * pow(x, 2);
}

double f2(double x){
    return std::sin(x);
}

double df2(double x){
    return std::cos(x);
}


TEST(NewtonRaphson, Comparison_equal)
{

double real1 = pow(2., 2./3);

double aux1 = NewtonRaphson(f1, df1, -5, 15, 5);
double aux2 = NewtonRaphson(f2, df2, -3./4 * M_PI, 1./2 * M_PI, -1./8 * M_PI);

ASSERT_NEAR(aux1, real1, 1.e-6 * real1);
ASSERT_NEAR(aux2, 0, 1.e-6 * aux2);

aux1 = NewtonRaphson(f1, df1, -5, 15, 15);
aux2 = NewtonRaphson(f2, df2, -3./4 * M_PI, 1./2 * M_PI, 1./2 * M_PI);

ASSERT_NEAR(aux1, real1, 1.e-6 * real1);
ASSERT_NEAR(aux2, 0, 1.e-6 * aux2);

}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

