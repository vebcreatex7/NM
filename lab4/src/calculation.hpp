#pragma once
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <iomanip>
#include <limits>
#define ld long double


ld eps(ld y1, ld y2);
ld Runge(ld y_half, ld y, int p);

//Task1
ld y_exact1(ld x);

//Euler
std::tuple<std::vector<ld>, std::vector<ld>> Euler_Method(ld h, ld x_0, ld x_1, ld y_0, ld z_0, std::ostream* log = nullptr);
std::vector<ld> Runge_Estimate_for_Euler(std::vector<ld> const & y, ld h, ld x_0, ld x_1, ld y_0, ld z_0);

//Runge-Kutta

std::tuple<std::vector<ld>, std::vector<ld>, std::vector<ld>> Runge_Kutta_Method(ld h, ld x_0, ld x_1, ld y_0, ld z_0, std::ostream* log = nullptr);
std::vector<ld> Runge_Estimate_for_Runge_Kutta(std::vector<ld> const & y, ld h, ld x_0, ld x_1, ld y_0, ld z_0);

//Adams
std::tuple<std::vector<ld>, std::vector<ld>> Adams_Method(ld h, ld x_0, ld x_1, ld y_0, ld z_0, std::ostream* log = nullptr);
std::vector<ld> Runge_Estimate_for_Adams(std::vector<ld> const & y, ld h, ld x_0, ld x_1, ld y_0, ld z_0);

//Task2

ld y_exact2(ld x);
std::tuple<std::vector<ld>, std::vector<ld>> Finite_Difference_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta, std::ostream* log = nullptr);
std::vector<ld> Runge_Estimate_for_Finite_Defference(std::vector<ld> const & y,ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta);

ld y_exact3(ld x);
std::tuple<std::vector<ld>, std::vector<ld>> Shooting_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta, std::ostream* log = nullptr);
std::vector<ld> Runge_Estimate_for_Shooting(std::vector<ld> const & y,ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta);