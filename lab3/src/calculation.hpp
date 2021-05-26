#pragma once
#include <string>
#include <iomanip>
#define ld long double
#include "matrix.hpp"
#include "Segment.hpp"


//Task1

std::tuple<std::string, ld, ld> Lagrange_Polynomial(std::vector<ld> const &, ld, std::ostream& log); 
std::tuple<std::string, ld, ld> Newton_Polynomial(std::vector<ld> const &, ld, std::ostream&);

//Task2

std::tuple<std::vector<TSegment>, ld> Cubic_Spline(std::vector<ld> const& x, std::vector<ld> const& f, ld Point);

//Task3
std::tuple<std::string, ld, std::string, ld> Least_Square_Method(std::vector<ld> const& x, std::vector<ld> const& y); 

//Task4

std::tuple<ld, ld> Numerical_Differentiation(std::vector<ld> const& x, std::vector<ld> const& y, ld Point);

//Task5

ld Rectangle_Method(ld x_0, ld x_k, ld h);
ld Trapezoid_Method(ld x_0, ld x_k, ld h);
ld Simpson_Method(ld x_0, ld x_k, ld h);
ld Runge_Romberg_Richardson_Method (ld F_half,ld F,ld p);