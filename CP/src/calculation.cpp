#include "calculation.hpp"


ld F(ld d2y, ld dy, ld y, ld x) {
    return d2y + log(dy) - std::exp(y) + x * x;
}
/*
ld A(std::vector<ld> const & y, ld x, size_t i, ld h) {
    ld part1 = F((y[i + 1] - 2. * y[i] + (y[i - 1] + delta)) / std::pow(h, 2.), (y[i + 1] - (y[i - 1] + delta)) / (2. * h), y[i], x);
    ld part2 = F((y[i + 1] - 2. * y[i] + y[i - 1]) / std::pow(h, 2.), (y[i + 1] - y[i - 1]) / (2. * h), y[i], x);
    return (part1 - part2) / delta;
}

ld B(std::vector<ld> const & y, ld x, size_t i, ld h) {
    ld part1 = F((y[i + 1] - 2. * (y[i] + delta) + y[i - 1]) / std::pow(h, 2.), (y[i + 1] - y[i - 1]) / (2. * h), (y[i] + delta), x);
    ld part2 = F((y[i + 1] - 2. * y[i] + y[i - 1]) / std::pow(h, 2.), (y[i + 1] - y[i - 1]) / (2. * h), y[i], x);
    return (part1 - part2) / delta;
}

ld C(std::vector<ld> const & y, ld x, size_t i, ld h) {
    ld part1 = F(((y[i + 1] + delta) - 2. * y[i] + y[i - 1]) / std::pow(h, 2.), ((y[i + 1] + delta) - y[i - 1]) / (2. * h), y[i], x);
    ld part2 = F((y[i + 1] - 2. * y[i] + y[i - 1]) / std::pow(h, 2.), (y[i + 1] - y[i - 1]) / (2. * h), y[i], x);
    return (part1 - part2) / delta;
}
*/




ld A(std::vector<ld> const & y, ld x, size_t i, ld h) {
    ld part1 = F((y[i + 1] - y[i]) / std::pow(h, 2.), (y[i + 1] - y[i]) / (2. * h), y[i], x);
    ld part2 = F((y[i + 1] - 2. * y[i] + y[i - 2]) / std::pow(h, 2.), (y[i + 1] - y[i - 2]) / (2. * h), y[i], x);
    return (part1 - part2) / (2. * h);
}

ld B(std::vector<ld> const & y, ld x, size_t i, ld h) {
    ld part1 = F((y[i - 1] - y[i + 1]) / std::pow(h, 2.), (y[i + 1] - y[i - 1]) / (2. * h), y[i + 1], x);
    ld part2 = F((y[i + 1] - y[i - 1]) / std::pow(h, 2.), (y[i + 1] - y[i - 1]) / (2. * h), y[i - 1], x);
    return (part1 - part2) / (2. * h);
}

ld C(std::vector<ld> const & y, ld x, size_t i, ld h) {
    ld part1 = F((y[i + 2] - 2. * y[i] + y[i - 1]) / std::pow(h, 2.), (y[i + 2] - y[i - 1]) / (2. * h), y[i], x);
    ld part2 = F((y[i - 1] - y[i]) / std::pow(h, 2.), (y[i] - y[i - 1]) / (2. * h), y[i], x);
    return (part1 - part2) / (2. * h);
}


TMatrix Jacobi_Matrix_Inverse(std::vector<ld> const & y, ld h, ld x_0, ld x_1) {
    size_t n = y.size() - 2;
    TMatrix J(n, n, 0.);
    J[0][0] = B(y, x_0 + h, 1, h);
    J[0][1] = C(y, x_0 + h, 1, h);

    for (size_t i = 1; i != n - 1; i++) {
        J[i][i - 1] = A(y, x_0 + (i + 1) * h, i + 1, h);
        J[i][i] = B(y, x_0 + (i + 1) * h, i + 1, h);
        J[i][i + 1] = C(y, x_0 + (i + 1) * h, i + 1, h);
    }
    J[n - 1][n - 2] = A(y, x_0 + n * h, n, h);
    J[n - 1][n - 1] = B(y, x_0 + n * h, n, h);

    auto [L, U, P] = J.LUdecomposition();

    return LU_Inverse_Matrix(L, U, P);
}

TMatrix F_column(std::vector<ld> const & y, ld h, ld x_0, ld x_1) {
    size_t n = y.size() - 2;
    TMatrix F_c(n, (size_t)1);
    for (size_t k = 0; k != n; k++) {
        size_t i = k + 1;
        F_c[k][0] = F((y[i + 1] - 2. * y[i] + y[i - 1]) / std::pow(h, 2), (y[i + 1] - y[i - 1]) / (2. * h), y[i], x_0 + i * h);
    }
    return F_c;
}



std::vector<ld> Newton_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta) {
    ld eps = 0.01;
    size_t n = (size_t)((x_1 - x_0) / h);

    std::vector<ld> y(n + 1);

    //initial approximation
    ld s = (x_1 + x_0) / (n + 1);
    y[1] = h;
    for (size_t i = 2; i != n; i++)
        y[i] = y[i - 1] + s;
    y[0] = alpha * h / (a_0 * h - a_1) - a_1 / (a_0 * h - a_1) * y[1];
    y[n] = beta * h  /(b_0 * h + b_1) + b_1 / (b_0 * h + b_1) * y[n - 1];


    TMatrix y_column(n - 1, (size_t)1);
        for (size_t i = 0; i != n - 1; i++)
            y_column[i][0] = y[i + 1];

    while (true) {

                
        TMatrix J_inv = Jacobi_Matrix_Inverse(y, h, x_0, x_1);
        TMatrix F_c = F_column(y, h, x_0, x_1);

        //Iteration
        TMatrix y_next = y_column - J_inv * F_c;


        
        for (size_t i = 0; i != n - 1; i ++)
            y[i + 1] = y_next[i][0];
        y[0] = alpha * h / (a_0 * h - a_1) - a_1 / (a_0 * h - a_1) * y[1];
        y[n] = beta * h  / (b_0 * h + b_1) + b_1 / (b_0 * h + b_1) * y[n - 1];

        if (TMatrix(y_next - y_column).Norm() < eps)
            break;
        
        y_column = y_next;
     
    }
    return y;

}

std::tuple<std::vector<ld>, std::vector<ld>> Finite_Difference_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta, std::ostream* log) {
    size_t n = (size_t)((x_1 - x_0) / h);
    std::vector<ld> x(n + 1);
    x[0] = x_0;
    for (size_t i = 1; i != n + 1; i++) 
        x[i] = x[i - 1] + h;
    std::vector<ld> y = Newton_Method(h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta);

    
    if (log) {
        for (auto a : x)
            *log << a << ' ';
        *log << '\n';
        for (auto a : y)
            *log << a << ' ';
        
    }

    std::cout << "Finite Difference Method is over\n";
    return std::make_tuple(x, y);
}



ld f1(ld x, ld y, ld z) {
    return z;
}

ld f2(ld x, ld y, ld z) {
    return -log(z) + std::exp(y) - x * x;
}


std::tuple<std::vector<ld>, std::vector<ld>, std::vector<ld>> Runge_Kutta_Method(ld h, ld x_0, ld x_1, ld y_0, ld z_0, std::ostream* log = nullptr) {
    size_t n = (size_t)((x_1 - x_0) / h);
    std::vector<ld> x(n + 1);
    std::vector<ld> y(n + 1);
    std::vector<ld> z(n + 1);

    x[0] = x_0;
    y[0] = y_0;
    z[0] = z_0;
    

    for (size_t i = 0; i != n; i++) {
        ld k_1 = h * f1(x[i], y[i], z[i]);
        ld l_1 = h * f2(x[i], y[i], z[i]);
        ld k_2 = h * f1(x[i] + h / 2., y[i] + k_1 / 2., z[i] + l_1 / 2.);
        ld l_2 = h * f2(x[i] + h / 2., y[i] + k_1 / 2., z[i] + l_1 / 2.);
        ld k_3 = h * f1(x[i] + h / 2., y[i] + k_2 / 2., z[i] + l_2 / 2.);
        ld l_3 = h * f2(x[i] + h / 2., y[i] + k_2 / 2., z[i] + l_2 / 2.);
        ld k_4 = h * f1(x[i] + h, y[i] + k_3, z[i] + l_3);
        ld l_4 = h * f2(x[i] + h, y[i] + k_3, z[i] + l_3);

        ld delta_y = (k_1 + (2. * k_2) + (2. * k_3) + k_4) / 6.;
        ld delta_z = (l_1 + (2. * l_2) + (2. * l_3) + l_4) / 6.;

        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + delta_y;
        z[i + 1] = z[i] + delta_z;

    }

    if (log) {

        for (auto a : x)
            *log << a << ' ';
        *log << '\n';

        for (auto a : y)
            *log << a << ' ';
        *log << '\n';
    }

    return std::make_tuple(x, y, z);
}


ld Phi (ld b_0, ld b_1, ld y_s, ld z_s, ld beta) {
    return b_0 * y_s + b_1 * z_s - beta;
}


std::tuple<std::vector<ld>, std::vector<ld>> Shooting_Method(ld h, ld x_0, ld x_1, ld a_0, ld a_1, ld alpha, ld b_0, ld b_1, ld beta, std::ostream* log) {
    std::vector<ld> x;
    std::vector<ld> y;

    ld eps = 0.01;
    int n = (int)((x_1 - x_0) / h);
    std::vector<ld> s;
    ld y_0, z_0;

    ld C_0, C_1;
    if (a_0 == 0.) {
        C_0 = -1. / a_1;
        C_1 = 0.;
    } else {
        C_0 = 0.;
        C_1 = -1. / a_0;
    }

    s.push_back((x_1 - x_0) / 2.);
    s.push_back(s[0] / 2.);

    std::vector<ld> y_s;
    std::vector<ld> z_s;

    
    y_0 = a_1 * s[0] - C_1 * alpha;
    z_0 = a_0 * s[0] - C_0 * alpha;
    auto ans1 = Runge_Kutta_Method(h, x_0, x_1, y_0, z_0);
        y_s.push_back(std::get<1>(ans1)[n]);
        z_s.push_back(std::get<2>(ans1)[n]);
    
    y_0 = a_1 * s[1] - C_1 * alpha;
    z_0 = a_0 * s[1] - C_0 * alpha;
    auto ans2 = Runge_Kutta_Method(h, x_0, x_1, y_0, z_0);
        y_s.push_back(std::get<1>(ans2)[n]);
        z_s.push_back(std::get<2>(ans2)[n]);



    

    size_t i = 2;
    while (true) {
        ld current_s = s[i - 1] - (s[i - 1] - s[i - 2]) / (Phi(b_0, b_1, y_s[i - 1], z_s[i - 1], beta) - Phi(b_0, b_1, y_s[i - 2], z_s[i - 2], beta)) * Phi(b_0, b_1, y_s[i - 1], z_s[i - 1], beta);
        s.push_back(current_s);
        y_0 = a_1 * s[i] - C_1 * alpha;
        z_0 = a_0 * s[i] - C_0 * alpha;
        auto ans = Runge_Kutta_Method(h, x_0, x_1, y_0, z_0);
        y_s.push_back(std::get<1>(ans)[n]);
        z_s.push_back(std::get<2>(ans)[n]);

        if (std::abs(Phi(b_0, b_1, y_s[i], z_s[i], beta)) < eps) {
            x = std::get<0>(ans);
            y = std::get<1>(ans);
            break;
        }

        i++;
    }

    if (log) {
        for (auto a : x)
            *log << a << ' ';
        *log << '\n';
        for (auto a : y)
            *log << a << ' ';
        
    }
    std::cout << "Shooting Method is over\n";
    return std::make_tuple(x, y);
}
