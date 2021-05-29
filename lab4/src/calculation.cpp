#include "calculation.hpp"
#include "matrix.hpp"


ld eps(ld y1, ld y2) {
    return std::abs(y1 - y2);
}


//Task1
ld Runge(ld y_half, ld y, int p) {
    return (y_half - y) / (std::pow(2., p) - 1.);
}

ld y_exact1(ld x) {
    return std::exp(x * x) + std::exp(x * std::sqrt(2)) + std::exp(-x * std::sqrt(2));
    
}


ld f1(ld x, ld y, ld z) {
    return z;
}

ld f2(ld x, ld y, ld z) {
    return 2. * y + 4. * x * x * std::exp(x * x);

}



//Euler

std::tuple<std::vector<ld>, std::vector<ld>> Euler_Method(ld h,  ld x_0, ld x_k, std::ostream* log) {
    size_t n = (size_t)((x_k - x_0) / h);
    std::vector<ld> x(n + 1);
    std::vector<ld> y(n + 1);
    std::vector<ld> z(n + 1);

    
    x[0] = 0.;
    y[0] = 3.;
    z[0] = 0.;
    

    for (size_t i = 0; i != n; i++) {
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + h * f1(x[i], y[i], z[i]);
        z[i + 1] = z[i] + h * f2(x[i], y[i], z[i]);
    }
    
    if (log) {
        for (auto a : x)
            *log << a << ' ';
        *log << '\n';

        for (auto a : y)
            *log << a << ' ';
        *log << '\n';
    }

    return std::make_tuple(x, y);
}


std::vector<ld> Runge_Estimate_for_Euler(std::vector<ld> const & y, ld x_0, ld x_k, ld h){
    auto [x_half, y_half] = Euler_Method(h / 2., x_0, x_k);
    size_t n = y.size();
    std::vector<ld> phi(n);

    for (size_t i = 0; i != n; i++) 
        phi[i] = Runge(y_half[2 * i], y[i], 1);
    
    return phi;

}


//Runge-Kutta


std::tuple<std::vector<ld>, std::vector<ld>, std::vector<ld>> Runge_Kutta_Method(ld h, ld x_0, ld x_k, std::ostream* log) {
    size_t n = (size_t)((x_k - x_0) / h);
    std::vector<ld> x(n + 1);
    std::vector<ld> y(n + 1);
    std::vector<ld> z(n + 1);
    
    x[0] = 0.;
    y[0] = 3.;
    z[0] = 0.;
    

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


std::vector<ld> Runge_Estimate_for_Runge_Kutta(std::vector<ld> const & y, ld x_0, ld x_k, ld h) {
    auto [x_half, y_half, other] = Runge_Kutta_Method(h / 2., x_0, x_k);
    size_t n = y.size();
    std::vector<ld> phi(n);

    for (size_t i = 0; i != n; i++)
        phi[i] = Runge(y_half[i], y[i], 4);

    return phi;
}


//Adams

std::tuple<std::vector<ld>, std::vector<ld>> Adams_Method(ld h, ld x_0, ld x_k, std::ostream* log) {
    size_t n = (size_t)((x_k - x_0) / h);
    auto [x_Runge, y_Runge, z_Runge] = Runge_Kutta_Method(h, x_0, 3. * h + x_0);

    std::vector<ld> f1_(n);
    std::vector<ld> f2_(n);
    std::vector<ld> x(n + 1);
    std::vector<ld> y(n + 1);
    std::vector<ld> z(n + 1);


    for (size_t i = 0; i != 4; i ++) {
        x[i] = x_Runge[i];
        y[i] = y_Runge[i];
        z[i] = z_Runge[i];
        f1_[i] = f1(x[i], y[i], z[i]);
        f2_[i] = f2(x[i], y[i], z[i]);
        
    }


    for (size_t i = 3; i != n; i++) {
        f1_[i] = f1(x[i], y[i], z[i]);
        f2_[i] = f2(x[i], y[i], z[i]);
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + h * (55. * f1_[i] - 59. * f1_[i - 1] + 37. * f1_[i - 2] - 9. * f1_[i - 3]) / 24.;
        z[i + 1] = z[i] + h * (55. * f2_[i] - 59. * f2_[i - 1] + 37. * f2_[i - 2] - 9. * f2_[i - 3]) / 24.;

    }

    if (log) {
        for (auto a : x)
            *log << a << ' ';
        *log << '\n';

        for (auto a : y)
            *log << a << ' ';
        *log << '\n';
    }



    return std::make_tuple(x, y);

}

std::vector<ld> Runge_Estimate_for_Adams(std::vector<ld> const & y, ld x_0, ld x_k, ld h) {
    auto [x_half, y_half] = Adams_Method(h / 2., x_0, x_k);
    size_t n = y.size();
    std::vector<ld> phi(n);

    for (size_t i = 0; i != n; i++)
        phi[i] = Runge(y_half[i], y[i], 4);

    return phi;

}





//Task2

ld y_exact2(ld x) {
    return 1. / x + 1.;
}

TMatrix Tridiagonal_Algorithm(TMatrix const& A, TMatrix const& D) {
    size_t n = A.Get_Rows();
    std::vector<long double> P(n);
    std::vector<long double> Q(n);

    //Forward Substitution

    long double a = 0.;
    long double b = A[0][1];
    long double c = A[0][2];
    long double d = D[0][0];
    P[0] = -c / b;
    Q[0] = d / b;
    for (size_t i = 1; i != n - 1; i++) {
        a = A[i][0];
        b = A[i][1];
        c = A[i][2];
        d = D[i][0];
        P[i] = -c / (b + a * P[i - 1]);
        Q[i] = (d - a * Q[i - 1]) / (b + a * P[i - 1]);
    }
    P[n - 1] = 0;
    a = A[n - 1][0];
    b = A[n - 1][1];
    d = D[n - 1][0];
    Q[n - 1] = (d - a * Q[n - 2]) / (b + a * P[n - 2]);

    

    //Back Substitution

    TMatrix x(n, (size_t)1);
    x[n - 1][0] = Q[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i][0] = P[i] * x[i + 1][0] + Q[i];

    return x;

}

ld p(ld x) {
    return 0.;
}

ld q(ld x) {
    return -2. / (x * x * (x + 1.));
}

ld f(ld x) {
    return 0.;
}




std::tuple<std::vector<ld>, std::vector<ld>> Finite_Difference_Method(ld h, ld x_0, ld x_k, std::ostream* log) {
    size_t n = (size_t)((x_k - x_0) / h);
    ld alpha = 0.;
    ld betta = -0.5;
    ld y_a = -1.;
    ld y_b = -1.;

    std::vector<ld> x(n + 1);
    x[0] = x_0;
    for (size_t i = 0; i != n; i++)
        x[i + 1] = x[i] + h;

    TMatrix A(n + 1, (size_t)3);
    TMatrix d(n + 1, (size_t)1);

    A[0][0] = 0.;
    A[0][1] = -2. / (h * (2 - p(x_0) * h)) + q(x_0) * h / (2. - p(x_0) * h) + alpha;
    A[0][2] = 2. / (h * (2 - p(x_0) * h));
    d[0][0] = y_a + h * f(x_0) / (2. - p(x_0) * h);

    for (size_t i = 1; i != n; i++) {
        A[i][0] = 1. / (h * h) - p(x[i]) / (2. * h);
        A[i][1] = -2. / (h * h) + q(x[i]);
        A[i][2] = 1. / (h * h) + p(x[i]) / (2. * h);
        d[i][0] = f(x[i]);
    }

    A[n][0] = -2. / (h * (2 + p(x[n]) * h));
    A[n][1] = 2 / (h * (2. + p(x[n]) * h)) - q(x[n]) * h / (2. + p(x[n]) * h) + betta;
    A[n][2] = 0.;
    d[n][0] = y_b - h * f(x[n]) / (2. + p(x[n]) * h);

    TMatrix Matrix_y = Tridiagonal_Algorithm(A, d);

    std::vector<ld> y(n + 1);
    for (size_t i = 0; i != n + 1; i++)
        y[i] = Matrix_y[i][0];

    if (log) {
        for (auto a : x)
            *log << a << ' ';
        *log << '\n';
        for (auto a : y)
            *log << a << ' ';
        
    }

    return std::make_tuple(x, y);


}
std::vector<ld> Runge_Estimate_for_Finite_Defference(std::vector<ld> const & y, ld x_0, ld x_k, ld h) {
    auto [x_half, y_half] = Finite_Difference_Method(h / 2., x_0, x_k);
    size_t n = y.size();
    std::vector<ld> phi(n);

    for (size_t i = 0; i != n; i ++)
        phi[i] = Runge(y_half[i], y[i], 1);
    
    return phi;
}
