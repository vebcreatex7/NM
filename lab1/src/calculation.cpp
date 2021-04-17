#include "calculation.hpp"
#include <vector>
#include <algorithm>


/*int Diagonals_Without_Zeros(TMatrix& A, TMatrix& b) {
    int p = 0;
    size_t n = A.Size();
    for (size_t i = 0; i != n; i++) {
        if (A[i][i] == 0) {
            for (size_t j = i + 1; j < n; j++) {
                if (A[j][i] != 0) {
                    A.Swap_Rows(i, j);
                    b.Swap_Rows(i, j);
                    p++;
                }
            }
        }
    }
    return p;
}*/

//Task1
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


long double LU_Determinant(TMatrix const& U, int p) {
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



//Task2
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



//Task3

std::tuple<TMatrix, int, int, long double> Iterative_Jacobi_Method(TMatrix const& A, TMatrix const& b, long double eps, std::ostream& log) {
    size_t n = A.Size();
    TMatrix alpha = A;
    TMatrix beta = b;
    for (size_t i = 0; i != n; i++) {
        if(alpha[i][i] == 0) {
            auto t = alpha.Change_Without_Zero(i);
            beta.Swap_Rows(t.first, t.second);
        }
        long double tmp = alpha[i][i];
        beta[i][0] /= tmp;
        for (size_t j = 0; j != n; j++) {
            
            if (j == i)
                alpha[i][i] = 0;
            else
                alpha[i][j] /= -tmp;

        }
    }
    //Norm
    long double norm = alpha.Norm();

    //A priori estimation of the number of iterations
    int k = (int)ceil(((log10(eps) - log10(beta.Norm()) + log10(1 - norm)) / log10(norm)) - 1);


    //Iterations

    TMatrix x_prev = beta;
    TMatrix x;
    int count = 0;
    log << "x_" << count << '\n' <<  x_prev << '\n';
    while (true) {
        count++;
        x = beta + alpha * x_prev;
        long double eps_i = TMatrix(x - x_prev).Norm();

        log << "x_" << count << '\n' <<  x_prev 
        << "eps_" << count << " = " << eps_i << "\n\n";

        if (eps_i <= eps)
            break;
        x_prev = x;
    }
    
    
    return std::make_tuple(x, count, k, norm);
}

std::tuple<TMatrix, int, int, long double> Seidel_Method(TMatrix const& A, TMatrix const& b, long double const eps, std::ostream& log) {
    size_t n = A.Size();
    TMatrix alpha = A;
    TMatrix beta = b;
    for (size_t i = 0; i != n; i++) {
        if(alpha[i][i] == 0) {
            auto t = alpha.Change_Without_Zero(i);
            beta.Swap_Rows(t.first, t.second);
        }
        long double tmp = alpha[i][i];
        beta[i][0] /= tmp;
        for (size_t j = 0; j != n; j++) {
            
            if (j == i)
                alpha[i][i] = 0;
            else
                alpha[i][j] /= -tmp;

        }
    }

    TMatrix B(n), C(n);
    for (size_t i = 0; i != n ; i++) {
        for (size_t j = 0; j != n; j++) {
            if (i > j)
                B[i][j] = alpha[i][j];
            else {
                C[i][j] = alpha[i][j];
            }
        }
    }
    TMatrix E(n);
    for (size_t i = 0; i != n; i++)
        E[i][i] = 1;

    alpha = (E - B).Inverse() * C;
    beta = (E - B).Inverse() * beta;

    
    //Norm
    long double norm = alpha.Norm();

    //A priori estimation of the number of iterations
    int k = (int)ceil(((log10(eps) - log10(beta.Norm()) + log10(1 - norm)) / log10(norm)) - 1);


    //Iterations

    TMatrix x_prev = beta;
    TMatrix x;
    int count = 0;
    log << "x_" << count << '\n' <<  x_prev << '\n';
    while (true) {
        count++;
        x = beta + alpha * x_prev;
        long double eps_i = TMatrix(x - x_prev).Norm();

        log << "x_" << count << '\n' <<  x_prev 
        << "eps_" << count << " = " << eps_i << "\n\n";

        if (eps_i <= eps)
            break;
        x_prev = x;
    }

    return std::make_tuple(x, count, k, norm);

}