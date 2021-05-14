#include "matrix.hpp"
#include "calculation.hpp"




void Task1() {
    std::fstream myFile;
    myFile.open("../tasks/t1.txt", std::ios::in);
    size_t n;
    myFile >> n;
    TMatrix A(n);
    TMatrix b(n, (size_t)1);
    myFile >> A >> b;
    myFile.close();

    auto LU = A.LUdecomposition();
    TMatrix L = std::get<0>(LU);
    TMatrix U = std::get<1>(LU);
    std::vector<std::pair<size_t, size_t>> P = std::get<2>(LU);
    TMatrix Origin = L * U;
    for (auto a : P) {
        if (a.first != a.second)
            Origin.Swap_Rows(a.first, a.second);
    }
    TMatrix x = LU_Solving_System(L, U, b, P);
    double det = LU_Determinant(U, P);
    TMatrix Inverse = LU_Inverse_Matrix(L, U, P);

    myFile.open("../ans/a1.txt", std::ios::out);
    myFile << "L:\n" << L
            << "\nU:\n" << U
            << "\nL*U:\n" << Origin
            << "\nx:\n" << x
            << "\n Det(A) = " << det
            << "\n\n A^(-1):\n" << Inverse
            << "\n A * A^(-1):\n" << A * Inverse;
    myFile.close();
    
}


void Task2() {
    std::fstream myFile;
    myFile.open("../tasks/t2.txt", std::ios::in);
    size_t n;
    myFile >> n;
    TMatrix A(n, (size_t)3);
    TMatrix d(n, (size_t)1);
    for (size_t i = 0; i != n; i++) {
        for (size_t j = 0; j != 3; j++) {
            if ((i == 0 and j == 0) or (i == n - 1 and j == 2))
                A[i][j] = 0;
            else {
                myFile >> A[i][j];
            }
        }
    }
    myFile >> d;
    myFile.close();

    TMatrix x = Tridiagonal_Algorithm(A, d);

    myFile.open("../ans/a2.txt", std::ios::out);
    myFile << "x:\n" << x;
    myFile.close();

}

void Task3() {
    
    std::fstream myFile;
    myFile.open("../tasks/t3.txt", std::ios::in);
    size_t n;
    long double eps;
    myFile  >> n;
    TMatrix A(n);
    TMatrix b(n, (size_t)1);
    myFile >> A >> b;
    myFile >> eps;
    myFile.close();

    std::fstream log;
    log.open("../logs/log3.1.txt", std::ios::out);
    auto ans = Iterative_Jacobi_Method(A, b, eps, log);
    log.close();

    TMatrix x = std::get<0>(ans);
    int count = std::get<1>(ans);
    int k = std::get<2>(ans);
    long double norm = std::get<3>(ans);

    
    myFile.open("../ans/a3.txt", std::ios::out);
    myFile  << "JACOBI METHOD\n\n"
            << "Norm of alpha:\n" << norm
            << "\nA priori estimation:\n" << k
            << "\nNumber of Iteration:\n" << count
            << "\nx:\n" << x;
    


    log.open("../logs/log3.2.txt", std::ios::out);
    ans = Seidel_Method(A, b, eps, log);
    log.close();


    x = std::get<0>(ans);
    count = std::get<1>(ans);
    k = std::get<2>(ans);
    norm = std::get<3>(ans);

    myFile  << "\n\nSEIDEL METHOD\n\n"
            << "Norm of alpha:\n" << norm
            << "\nA priori estimation:\n" << k
            << "\nNumber of Iteration:\n" << count
            << "\nx:\n" << x;
    myFile.close();
    
}


void Task4() {
    std::fstream myFile;
    myFile.open("../tasks/t4.txt", std::ios::in);
    size_t n;
    long double eps;
    myFile >> n;
    TMatrix A(n);
    myFile >> A;
    myFile  >> eps;
    myFile.close();

    std::fstream log;
    log.open("../logs/log4.txt", std::ios::out);
    auto ans = Rotation_Method(A, eps, log);
    log.close();

    
    myFile.open("../ans/a4.txt", std::ios::out);
    std::vector<long double> Eigenvalues = std::get<0>(ans);
    std::vector<TMatrix> Eigenvectors = std::get<1>(ans);
    size_t count = std::get<2>(ans);
    for (size_t i = 0; i != n; i++) {
        myFile << "Eigenvalue_" << i + 1 << " = " << Eigenvalues[i] << '\n'
                << "x_" << i + 1 << ":\n" << Eigenvectors[i]
                << "A * x - labda * x: \n" << A * Eigenvectors[i] - Eigenvectors[i] * Eigenvalues[i] << "\n\n";
    }
    myFile << "Number of iterations: " << count << '\n';
    myFile.close();
}

void Test() {
    size_t n;
    std::cin >> n;
    TMatrix t(n);
    std::cin >> t;

    auto [L, U, P] = t.LUdecomposition();
    std::cout << L * U;
    for (auto a : P)
        std::cout << a.first << ' ' << a.second << '\n';
}



int main() {
    
    Task1();
    Task2();
    Task3();
    Task4();
    //Test();
    return 0;

}