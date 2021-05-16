#include "calculation.hpp"



//Task1
long double f(long double x) {
    return sqrt(1 - std::pow(x, 2.)) - exp(x) + 0.1;
}

long double First_Derivative(long double x) {
    return -x / sqrt(1 - std::pow(x, 2.)) - exp(x);
}

/*
long double Second_Derivative(long double x) {
    return -std::pow(x, 2.) / std::pow((1 - x * x), 3/2) - 1 / sqrt(1 - x * x) - exp(x);
}
*/


std::tuple<long double, int> Newton_Method(long double a, long double b, long double eps, std::ostream& log) {
    long double x_prev = b;
    long double x;
    int count = 0;
    do {
        count++;
        x = x_prev - f(x_prev) / First_Derivative(x_prev);
        log << std::setprecision(7) << std::fixed  << x << ' ' << f(x) << ' ' << First_Derivative(x) << ' ' << - f(x) / First_Derivative(x) << '\n';

    } while(std::abs(x - x_prev) >= eps && (x_prev = x));

    return std::make_tuple(x, count);

}


long double phi(long double x) {
    return log(0.1 + std::sqrt(1 - x * x));
}

/*
long double First_Derivative_Phi(long double x) {
    return x / (x * x - 0.1 * std::sqrt(1 - x * x) - 1);
}
*/


std::tuple<long double, int> Simple_Iterations_Method(long double a, long double b, long double eps, std::ostream& log) {
    long double q = 0.6; // |phi'(x) < 0.6|
    long double x_prev = (a + b) / 2.;
    long double x;
    int count = 0;
    do
    {
        count++;
        x = phi(x_prev);
        log << std::setprecision(7) << std::fixed  << x << ' ' << f(x) << ' ' << First_Derivative(x) << ' ' << - f(x) / First_Derivative(x) << '\n';
        
    } while (q/(1 - q) * std::abs(x - x_prev) >= eps && (x_prev = x));
    
    return std::make_tuple(x, count);
}