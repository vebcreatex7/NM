#pragma once
#include <iostream>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

#define ld long double

struct TSegment {
    ld a, b, c, d;
    void Print(std::ostream& os) {
        os << std::setprecision(5) << std::fixed << a << ' ' << b << ' ' << c << ' ' << d << '\n';
    };
};


