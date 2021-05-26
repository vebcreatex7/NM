#include "calculation.hpp"
#include <algorithm>




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




//Task1

ld f(ld x) {
    return tan(x);
}


ld omega_t(size_t pivot, std::vector<ld> const & v) {
    ld val = 1.;
    for (size_t i = 0; i != v.size(); i++) {
        if (i != pivot)
            val *= (v[pivot] - v[i]);
    }
    return val;
}


std::tuple<std::string, ld, ld> Lagrange_Polynomial(std::vector<ld>  const & points, ld Exact_X,  std::ostream& log) {
    ld value = 0.;

    std::string polynomial = "";
    for (size_t i = 0; i != points.size(); i++) {
        ld f_i = f(points[i]);
        ld omega_i = omega_t(i, points);
        log << std::setprecision(5) << std::fixed << i << ' ' << points[i] << ' ' << f_i << ' ' << omega_i << ' ' << f_i / omega_i << ' ' << Exact_X - points[i] << '\n';
        if (f_i / omega_i == 0)
            continue;
        if (f_i / omega_i > 0 && i > 0)
            polynomial += " + ";
        polynomial += std::to_string(f_i / omega_i);
        ld tmp = 1.;
        for (size_t j = 0; j != points.size(); j++) {
            if (j != i) {
                if (points[j] == 0)
                    polynomial += "x";
                else
                    polynomial += "(x - " + std::to_string(points[j]) + ")";
                tmp *= (Exact_X - points[j]);
            }
        }

        value += f_i / omega_i *  tmp;
    }


    return std::make_tuple(polynomial, value, f(Exact_X));
}


std::vector<std::vector<ld>> Divided_Difference(std::vector<ld> const & points) {
    size_t n = points.size();
    std::vector<std::vector<ld>> f_i(n, std::vector<ld>(n, 0));
    for (size_t i = 0; i != n; i++)
        f_i[0][i] = f(points[i]);
    for (size_t i = 1; i != n; i++) {
        for (size_t j = 0; j != n - i; j++)
            f_i[i][j] = (f_i[i - 1][j] - f_i[i - 1][j + 1]) / (points[j] - points[j + i]);
    }

    return f_i;
}



std::tuple<std::string, ld, ld> Newton_Polynomial(std::vector<ld> const & points, ld Exact_X, std::ostream& log) {
    size_t n = points.size();
    std::string polynomial = "";
    ld value = 0.;
    std::vector<std::vector<ld>> f_i = Divided_Difference(points);
    
    for (size_t i = 0; i != n; i++) {
        log << std::setprecision(5) << std::fixed << i << ' ' << points[i] << ' ';
        for (size_t j = 0; j != n - i; j++)
            log << std::setprecision(5) << std::fixed << f_i[j][i] << ' ';
        log << '\n';
    }
    
    for (size_t i = 0; i != n; i++) {
        ld tmp = 1.;
        if (f_i[i][0] == 0)
            continue;
        else {
            if (f_i[i][0] > 0 && i > 0)
                polynomial += " + "; 
            polynomial += std::to_string(f_i[i][0]);
            for (size_t j = 0; j != i; j ++) {
                tmp *= (Exact_X - points[j]);
                if (points[j] == 0)
                    polynomial += "x";
                else
                    polynomial += "(x - " + std::to_string(points[j]) + ")";
            }
            value += f_i[i][0] * tmp;
            
        }
    }
    return std::make_tuple(polynomial, value, f(Exact_X));
}


//Task2



ld h(std::vector<ld> const& x, size_t i) {
    return x[i] - x[i - 1];
}

std::tuple<std::vector<TSegment>, ld> Cubic_Spline(std::vector<ld> const& x, std::vector<ld> const& f, ld Point) {
    size_t n = x.size();
    n -= 2;
    TMatrix A(n, (size_t)3);
    TMatrix d(n, (size_t)1);

    A[0][0] = 0;
    A[0][1] = 2. * (h(x, 1) + h(x, 2));
    A[0][2] = h(x, 2);
    for (size_t i = 1; i < n - 1; ++i) {
        A[i][0] = h(x, i+1);
        A[i][1] = 2.0 * (h(x, i+1) + h(x, i+2));
        A[i][2] = h(x, i+2);
    }
    A[n - 1][0] = h(x, n - 2);
    A[n - 1][1] = 2. * (h(x, n - 2) + h(x, n - 1));
    A[n - 1][2] = 0;

    for (size_t i = 0; i != n; i++)
        d[i][0] = 3.0 * ((f[i+2] - f[i+1])/h(x, i+2) - (f[i+1] - f[i]) / h(x, i+1));
    
    

    TMatrix cc = Tridiagonal_Algorithm(A, d);
    std::vector<ld> c(n + 1);
    c[0] = 0.;
    for (size_t i = 0; i != n; i++)
        c[i + 1] = cc[i][0];

    n += 2;

    std::vector<TSegment> seg_(n - 1);
    for (size_t i = 0; i < n - 2; ++i) {
        seg_[i].a = f[i];
        seg_[i].b = (f[i+1] - f[i]) / h(x, i+1) - h(x, i+1)*(c[i+1] + 2.0 * c[i]) / 3.0;
        seg_[i].c = c[i];
        seg_[i].d = (c[i+1] - c[i]) / (3.0 * h(x, i+1));
    }

    seg_[n-2].a = f[n-2];
    seg_[n-2].b = (f[n-1] - f[n-2]) / h(x, n-1) - 2.0 * h(x, n-1) * c[n-2] / 3.0;
    seg_[n-2].c = c[n-2];
    seg_[n-2].d = -c[n-2]/ (3.0 * h(x, n-1));
    

    auto it = std::lower_bound(x.begin(), x.end(), Point);
    size_t  b = it - x.begin();
    size_t a = b - 1;
    TSegment s = seg_[a];
    ld h = (Point - x[a]);
    ld f_x = s.a + s.b * h + s.c * h * h + s.d * h * h * h;
    std::vector<std::vector<ld>> v;
    return std::make_tuple(seg_, f_x);

    
}


//Task4

ld First_Derivative(std::vector<ld> const& x, std::vector<ld> const& y, ld Point, size_t i) {
    ld part1 = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    ld part2 = ((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - part1) / (x[i + 2] - x[i]) * (2. * Point - x[i] - x[i + 1]);
    return part1 + part2;
}

ld Second_Derivative(std::vector<ld> const& x, std::vector<ld> const& y, ld Point, size_t i) {
    ld part1 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    ld part2 = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    return 2. * (part1 - part2) / (x[i + 2] - x[i]);
}

std::tuple<ld, ld> Numerical_Differentiation(std::vector<ld> const& x, std::vector<ld> const& y, ld Point) {
    auto it = std::lower_bound(x.begin(), x.end(), Point);
    size_t i = it - x.begin() - 1;
    return std::make_tuple(First_Derivative(x, y, Point, i), Second_Derivative(x, y, Point, i));



}


std::tuple<TMatrix, std::string> Approximating_Polynomial_First_Degree(std::vector<ld> const& x, std::vector<ld> const& y) {
    size_t d = 2;
    TMatrix A(d);
    TMatrix b(d, (size_t)1);
    size_t n = x.size();
    A[0][0] = n;
    for (auto a : x)
        A[0][1] += a;
    A[1][0] = A[0][1];
    for (auto a : x)
        A[1][1] += a * a;

    for (auto a : y)
        b[0][0] += a;
    for (size_t i = 0; i != x.size(); i++)
        b[1][0] += x[i] * y[i];


    auto [L, U, P] = A.LUdecomposition();
    
    TMatrix a = LU_Solving_System(L, U, b, P);

    std::string Polynomial = std::to_string(a[0][0]);
    if (a[1][0] > 0.)
        Polynomial += " + ";
    Polynomial += std::to_string(a[1][0]) + "x";

    return std::make_tuple(a, Polynomial);

}

std::tuple<TMatrix, std::string> Approximating_Polynomial_Second_Degree(std::vector<ld> const& x, std::vector<ld> const& y) {
    size_t d = 3;
    TMatrix A(d);
    TMatrix b(d, (size_t)1);

    A[0][0] = (ld)x.size();
    for (auto a : x) {
        A[0][1] += a;
        A[0][2] += a * a;
        A[1][2] += a * a * a;
        A[2][2] += a * a * a * a;
    }

    A[1][0] = A[0][1];

    A[1][1] = A[0][2];
    A[2][0] = A[1][1];

    A[2][1] = A[1][2];

    for (size_t i = 0; i != y.size(); i++) {
        b[0][0] += y[i];
        b[1][0] += x[i] * y[i];
        b[2][0] += x[i] * x[i] * y[i];
    }

    auto [L, U, P] = A.LUdecomposition();
    TMatrix a = LU_Solving_System(L, U, b, P);

    std::string Polynomial = std::to_string(a[0][0]);
    if (a[1][0] > 0.)
        Polynomial += " + ";
    Polynomial += std::to_string(a[1][0]) + "x";

    if (a[2][0] > 0.)
        Polynomial += " + ";
    Polynomial += std::to_string(a[2][0]) + "x^2";

    return std::make_tuple(a, Polynomial);
        
}

ld F_1(TMatrix a, ld x) {
    return a[0][0] + a[1][0] * x;
}

ld F_2(TMatrix a, ld x) {
    return a[0][0] + a[1][0] * x + a[2][0] * x * x;
}

std::tuple<std::string, ld, std::string, ld> Least_Square_Method(std::vector<ld> const& x, std::vector<ld> const& y) {
    size_t n = x.size();
    auto [a1, p1] = Approximating_Polynomial_First_Degree(x, y);

    ld PHI1 = 0.;
    for (size_t i = 0; i != n; i++)
        PHI1 += std::pow(F_1(a1, x[i]) - y[i], 2.);
    
    auto [a2, p2] = Approximating_Polynomial_Second_Degree(x, y);

    ld PHI2 = 0;
    for (size_t i = 0; i != n; i++)
        PHI2 += std::pow(F_2(a2, x[i]) - y[i], 2.);

    return std::make_tuple(p1, PHI1, p2, PHI2);

}


//Task5

ld y(ld x) {
    return x / std::pow((3. * x + 4.), 3.);
}

ld Rectangle_Method(ld x_0, ld x_k, ld h) {
    int n = int(x_k - x_0) / h;
    ld F = 0.;
    
    ld x_cur = x_0;
    for (int i = 0; i != n - 1; i++) {
        F += y((x_cur + x_cur + h) / 2.);
        x_cur += h;
    } 
    F *= h;
    return F;
}


ld Trapezoid_Method(ld x_0, ld x_k, ld h) {
    int n = int(x_k - x_0) / h;
    ld F = 0.;

    ld x_cur = x_0;
    for (int i = 0; i != n - 1; i++) {
        F += (y(x_cur) + y(x_cur + h));
        x_cur += h;
    }

    F *= (h / 2.);

    return F;
}


ld Simpson_Method(ld x_0, ld x_k, ld h) {
    int n = int(x_k - x_0) / h;
    ld F = 0.;

    ld x_cur = x_0;
    for (int i = 0; i != n - 2; i++) {
        F += (y(x_cur) + 4 * y(x_cur + h) + y(x_cur + 2 * h));
        x_cur +=2 * h;
    }
    F *= (h / 3.);
    return F;
}

ld Runge_Romberg_Richardson_Method (ld F_half,ld F,ld p) {
    return (F_half - F) / (std::pow(2, p) - 1);
}