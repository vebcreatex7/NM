#pragma once
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <iomanip>
#include <limits>
#include "matrix.hpp"
#define ld long double


std::tuple<std::vector<ld>, std::vector<ld>> Shooting_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta, std::ostream* log = nullptr);
std::tuple<std::vector<ld>, std::vector<ld>> Finite_Difference_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta, std::ostream* log = nullptr);

