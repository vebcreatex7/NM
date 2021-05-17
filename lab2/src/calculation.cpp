#include "calculation.hpp"
#include "matrix.hpp"



TMatrix LU_Solving_System(TMatrix const& L, TMatrix const& U, TMatrix b, std::vector<std::pair<size_t, size_t>> const& p) {
    for (size_t i = 0; i != p.size(); i++)
        b.Swap_Rows(p[i].first, p[i].second);
    //Ly = b
    //Forward Substitution

    size_t n = L.Size();
    TMatrix y(n, size_t(1));
    for (size_t i = 0; i != n; i++) {
        long double t = 0.;
        for (size_t j = 0; j != i; j++)
           t += y[j][0] * L[i][j];
        y[i][0] = b[i][0] - t;
    }
    
    //Ux = y;
    //Back Substitution

    TMatrix x(n, (size_t)1);
    for (int i = n - 1; i >= 0; i--) {
        long double t = 0.;
        for (int j = i + 1; j < (int)n; j++)
            t += U[i][j] * x[j][0];
        x[i][0] = (y[i][0] - t) / U[i][i];
    }

    return x;

}


long double LU_Determinant(TMatrix const& U, std::vector<std::pair<size_t, size_t>> const & P) {
    size_t p = 0;
    for (auto a : P)
        p = a.first != a.second ? p + 1 : p;

    long double det = 1.;
    for (size_t i = 0; i != U.Size(); i++)
        det *= U[i][i];
    return std::pow(-1, p) * det; 
}



TMatrix LU_Inverse_Matrix(TMatrix const& L, TMatrix const& U, std::vector<std::pair<size_t, size_t>> const& p) {
    size_t n = L.Size();
    TMatrix Inverse(n);
    for (size_t i = 0; i != n; i++) {
        TMatrix b(n, (size_t)1);
        b[i][0] = 1;
        TMatrix tmp = LU_Solving_System(L, U, b, p);
        for (size_t j = 0; j != n; j++)
            Inverse[j][i] = tmp[j][0];
    }
    return Inverse;
}

int Sign(long double d) {
    return d < 0 ? -1 : 1;
}






//Task1

//Newton
long double f(long double x) {
    return sqrt(1 - std::pow(x, 2.)) - exp(x) + 0.1;
}

long double First_Derivative(long double x) {
    return -x / sqrt(1 - std::pow(x, 2.)) - exp(x);
}

/*
long double Second_Derivative(long double x) {
    return -std::pow(x, 2.) / std::pow((1 - x * x), 3/2) - 1 / sqrt(1 - x * x) - exp(x);
}
*/



std::tuple<long double, int> Newton_Method(long double a, long double b, long double eps, std::ostream& log) {
    long double x_prev;
    long double x = b;
    int count = 0;
    do {
        x_prev = x;
        count++;
        x = x_prev - f(x_prev) / First_Derivative(x_prev);
        log << std::setprecision(7) << std::fixed  << x << ' ' << f(x) << ' ' << First_Derivative(x) << ' ' << - f(x) / First_Derivative(x) << '\n';

    } while(std::abs(x - x_prev) >= eps);

    return std::make_tuple(x, count);

}




//Iterations
long double phi(long double x) {
    return log(0.1 + std::sqrt(1 - x * x));
}

/*
long double First_Derivative_Phi(long double x) {
    return x / (x * x - 0.1 * std::sqrt(1 - x * x) - 1);
}
*/


std::tuple<long double, int> Simple_Iterations_Method(long double a, long double b, long double eps, std::ostream& log) {
    long double q = 0.6; // |phi'(x)| < 0.6
    long double x_prev;
    long double x  = (a + b) / 2.;
    int count = 0;
    do
    {
        x_prev = x;
        count++;
        x = phi(x_prev);
        log << std::setprecision(7) << std::fixed  << x << ' ' << f(x) << ' ' << First_Derivative(x) << ' ' << - f(x) / First_Derivative(x) << '\n';
        
    } while (q/(1 - q) * std::abs(x - x_prev) >= eps);
    
    return std::make_tuple(x, count);
}




//Task2


//Newton
long double f1(long double x1, long double x2) {
    return (x1 * x1 + 16.) * x2 - 64.;
}

long double f2(long double x1, long double x2) {
    return std::pow((x1 - 2.), 2.) + std::pow((x2 - 2.), 2.) - 16;
}

long double df1_dx1(long double x1, long double x2) {
    return 2 * x1 * x2;
}

long double df1_dx2(long double x1, long double x2) {
    return x1 * x1 + 16;
}

long double df2_dx1(long double x1, long double x2) {
    return 2 * (x1 - 2);
}
long double df2_dx2(long double x1, long double x2) {
    return 2 * (x2 - 2);
}

TMatrix Jacobi_Matrix(long double x1, long double x2) {
    TMatrix J(2);
    J[0][0] = df1_dx1(x1, x2);
    J[0][1] = df1_dx2(x1, x2);
    J[1][0] = df2_dx1(x1, x2);
    J[1][1] = df2_dx2(x1, x2);

    return J;
}

TMatrix A_1 (long double x1, long double x2) {
    TMatrix A1(2);
    A1[0][0] = f1(x1, x2);
    A1[0][1] = df1_dx2(x1, x2);
    A1[1][0] = f2(x1, x2);
    A1[1][1] = df2_dx2(x1, x2);
    return A1;
}

TMatrix A_2 (long double x1, long double x2) {
    TMatrix A2(2);
    A2[0][0] = df1_dx1(x1, x2);
    A2[0][1] = f1(x1, x2);
    A2[1][0] = df2_dx1(x1, x2);
    A2[1][1] = f2(x1, x2);
    return A2;
}

long double Square_Determinant(TMatrix const& M) {
    return M[0][0] * M[1][1] - M[0][1] * M[1][0];
}

std::tuple<TMatrix, int> Dimentional_Newton_Method(long double eps, std::ostream& log) {
    TMatrix x_prev(2, (size_t)1);
    TMatrix x(2, (size_t)1);
    x[0][0] = 5.5;
    x[1][0] = 2.;

    int count = 0;

    do {
        count++;
        x_prev = x;
        TMatrix J = Jacobi_Matrix(x_prev[0][0], x_prev[1][0]);
        TMatrix A1 = A_1(x_prev[0][0], x_prev[1][0]);
        TMatrix A2 = A_2(x_prev[0][0], x_prev[1][0]);

        long double det_J = Square_Determinant(J);
        long double det_A1 = Square_Determinant(A1);
        long double det_A2 = Square_Determinant(A2);

        x[0][0] = x_prev[0][0] - det_A1 / det_J;
        x[1][0] = x_prev[1][0] - det_A2 / det_J;

        log << "\nk = " << count << '\n' << "x_1 = " << x[0][0] << " x_2 = " << x[1][0] << '\n';

    } while(TMatrix(x - x_prev).Norm() >= eps);

    return std::make_tuple(x, count);

}






//Iterations

long double phi1(long double x1 = 0., long double x2 =  0.) {
    return std::sqrt(16 - std::pow((x2 - 2), 2.)) + 2;
}

long double phi2(long double x1, long double x2 = 0.) {
    return 64. / (x1 * x1 + 16);
}

long double dphi1_dx1 (long double x1 = 0., long double x2 =  0.) {
    return 0.;
}

long double dphi1_dx2 (long double x1 = 0., long double x2 =  0.) {
    return (x2 - 2.) / (std::sqrt(-x2 * x2 + 4 * x2 + 12));;
}

long double dphi2_dx1 (long double x1 = 0., long double x2 =  0.) {
    return -128 * x1 / std::pow((x1 * x1 + 16), 2);
}

long double dphi2_dx2 (long double x1 = 0., long double x2 =  0.) {
    return 0.;
}


std::tuple<TMatrix, int> Dimentional_Simple_Iterations_Method(long double eps, std::ostream& log) {
    long double q = 0.5;
    TMatrix x_prev(2, size_t(1));
    TMatrix x(2, (size_t)1);
    x[0][0] = 5.5;
    x[1][0] = 2.;

    int count = 0;
    do
    {
        x_prev = x;
        x[0][0] = phi1(x_prev[0][0], x_prev[1][0]);
        x[1][0] = phi2(x_prev[0][0], x_prev[1][0]);
        count++;

        log << "\nk = " << count << '\n' << "x_1 = " << x[0][0] << " x_2 = " << x[1][0] << '\n';

    } while (q / (1 - q) * TMatrix(x - x_prev).Norm() >= eps);

    return std::make_tuple(x, count);
    
}


