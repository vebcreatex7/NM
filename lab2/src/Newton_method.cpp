#include "Newton_method.hpp"

TNewton_Method::TNewton_Method(long double a, long double b)
    : a(a), b(b) {}

long double TNewton_Method::f(long double x) {
    return sqrt(1 - std::pow(x, 2.)) - exp(x) + 0.1;
}

long double TNewton_Method::First_Derivative(long double x) {
    return -x / sqrt(1 - std::pow(x, 2.)) - exp(x);
}

long double TNewton_Method::Second_Derivative(long double x) {
    return -std::pow(x, 2.) / std::pow((1 - x * x), 3/2) - 1 / sqrt(1 - x * x) - exp(x);
}

