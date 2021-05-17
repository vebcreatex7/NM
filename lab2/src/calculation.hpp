#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <iomanip>

#include "matrix.hpp"


TMatrix LU_Solving_System(TMatrix const& L, TMatrix const& U, TMatrix b, std::vector<std::pair<size_t, size_t>> const& p);
long double LU_Determinant(TMatrix const& U, std::vector<std::pair<size_t, size_t>> const&);
TMatrix LU_Inverse_Matrix(TMatrix const& L, TMatrix const& U, std::vector<std::pair<size_t, size_t>> const& p);
TMatrix Seidel_Method(TMatrix const& A, TMatrix const& b, long double const eps);
int Sign(long double d);



//Task1


std::tuple<long double, int> Newton_Method(long double a, long double b, long double eps, std::ostream& log);
std::tuple<long double, int> Simple_Iterations_Method(long double a, long double b, long double eps, std::ostream& log);


//Task2

std::tuple<TMatrix, int> Dimentional_Newton_Method(long double eps, std::ostream& log);
std::tuple<TMatrix, int> Dimentional_Simple_Iterations_Method(long double eps, std::ostream& log);

