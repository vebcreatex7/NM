#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <iomanip>



//Task1

long double f(long double x);
long double First_Derivative(long double x);
long double Second_Derivative(long double x);
std::tuple<long double, int> Newton_Method(long double a, long double b, long double eps, std::ostream& log);

long double phi(long double x);
long double First_Derivative_Phi(long double x);
std::tuple<long double, int> Simple_Iterations_Method(long double a, long double b, long double eps, std::ostream& log);