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