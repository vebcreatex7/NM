#include "calculation.hpp"
#include "matrix.hpp"


void task1() {

    //a
    size_t n = 4;
    std::vector<ld> v1(n);
    v1[0] = 0, v1[1] = M_PI / 8, v1[2] = 2 * M_PI / 8, v1[3] = 3 * M_PI / 8;

    std::vector <ld> v2(n);
    v2[0] = 0, v2[1] = M_PI / 8, v2[2] = M_PI / 3, v2[3] = 3 * M_PI / 8;

    ld Exact_X = 3 * M_PI / 16;
    
    std::fstream log1("../logs/log1.1.txt", std::ios::out);
    auto [Polynomial1, L_x, f1_x] = Lagrange_Polynomial(v1, Exact_X, log1);
    log1.close();

    std::fstream myFile("../ans/a1.1.txt", std::ios::out);
    myFile << std::setprecision(25) << "Lagrange Polynomial:\n" << "L(x) = " << Polynomial1 << "\n" 
        << "L(X*) = " << L_x << ' ' << "f(X*) = " << f1_x
        << "\nAbsolute error of interpolation:" << std::abs(L_x - f1_x) << "\n\n";
    
    log1.open("../logs/log1.2.txt");
    auto [Polynomial2, P_x, f_x] = Newton_Polynomial(v1, Exact_X, log1);
    log1.close();
    myFile << "Newton Polynomial:\n" << "P(x) = " << Polynomial2 << "\n" 
        << "P(X*) = " << P_x << ' ' << "f(X*) = " << f_x
        << "\nAbsolute error of interpolation:" << std::abs(P_x - f_x) << "\n\n";
    myFile.close();

    // b

    log1.open("../logs/log1.3.txt");
    auto [Polynomial3, L_x1, f1_x1] = Lagrange_Polynomial(v2, Exact_X, log1);
    log1.close();

    myFile.open("../ans/a1.2.txt");
    myFile << std::setprecision(25) << "Lagrange Polynomial:\n" << "L(x) = " << Polynomial3 << "\n" 
        << "L(X*) = " << L_x1 << ' ' << "f(X*) = " << f1_x1
        << "\nAbsolute error of interpolation:" << std::abs(L_x1 - f1_x1) << "\n\n";

    log1.open("../logs/log1.4.txt");
    auto [Polynomial4, P_x1, f_x1] = Newton_Polynomial(v2, Exact_X, log1);
    log1.close();

    myFile << "Newton Polynomial:\n" << "P(x) = " << Polynomial4 << "\n" 
        << "P(X*) = " << P_x1 << ' ' << "f(X*) = " << f_x1
        << "\nAbsolute error of interpolation:" << std::abs(P_x1 - f_x1) << "\n\n";
    myFile.close();
    


}


void task2() {
    std::fstream myFile("../tasks/t2.txt", std::ios::in);
    size_t n;
    myFile >> n;
    std::vector<ld> x(n);
    std::vector<ld> f(n);
    for (auto & a : x)
        myFile >> a;
    for (auto & a : f)
        myFile >> a;

    ld Point;
    myFile >> Point;
    myFile.close();

    auto [Segments, f_x] = Cubic_Spline(x, f, Point);
    
    std::fstream log("../logs/log2.txt", std::ios::out);
    for (auto a : Segments)
        a.Print(log);
    
    log.close();
    
    myFile.open("../ans/a2.txt", std::ios::out);
    myFile << "f(x) = " << f_x << '\n';


}


void task4() {
    std::fstream myFile("../tasks/t4.txt", std::ios::in);
    size_t n;
    myFile >> n;
    std::vector<ld> x(n);
    std::vector<ld> y(n);
    for (auto & a : x)
        myFile >> a;
    for (auto & a : y)
        myFile >> a;

    ld Point;
    myFile >> Point;
    myFile.close();


    myFile.open("../ans/a4.txt", std::ios::out);
    auto [d1, d2] =  Numerical_Differentiation(x, y, Point);
    myFile << d1 << ' ' << d2 << '\n';
    myFile.close();

}

void task3() {
    std::fstream myFile("../tasks/t3.txt", std::ios::in);
    size_t n;
    myFile >> n;
    std::vector<ld> x(n);
    std::vector<ld> y(n);
    for (auto & a : x)
        myFile >> a;
    for (auto & a : y)
        myFile >> a;

    myFile.close();


    auto [p1, PHI1, p2, PHI2] = Least_Square_Method(x, y);

    myFile.open("../ans/a3.txt", std::ios::out);

    myFile << "F1(x) = " << p1 << '\n'
            << "PHI = " << PHI1 << '\n'
            << "F2(x) = " << p2 << '\n'
            << "PHI = " << PHI2 << '\n';

    myFile.close();
}

void task5() {
    std::fstream myFile("../tasks/t5.txt", std::ios::in);
    ld x_0, x_k, h1, h2;
    myFile >> x_0 >> x_k >> h1 >> h2;
    myFile.close();

    ld Exact_Solution = -0.12245;

    ld R1 = Rectangle_Method(x_0, x_k, h1);
    ld T1 = Trapezoid_Method(x_0, x_k, h1);
    ld S1 = Simpson_Method(x_0, x_k, h1);

    ld R2 = Rectangle_Method(x_0, x_k, h2);
    ld T2 = Trapezoid_Method(x_0, x_k, h2);
    ld S2 = Simpson_Method(x_0, x_k, h2);

    ld posteriori_R = Runge_Romberg_Richardson_Method(R2, R1, 1);
    ld posteriori_T = Runge_Romberg_Richardson_Method(T2, T1, 2);
    ld posteriori_S = Runge_Romberg_Richardson_Method(S2, S1, 3);



    myFile.open("../ans/a5.txt", std::ios::out);
    myFile  << "Exact Solution: " << Exact_Solution << "\n\n"
            << "h1 = " << h1 << '\n'
            << "Rectangle Method: " << R1 <<'\n'
            << "Trapezoid Method: " << T1 << '\n'
            << "Simpson Method: " << S1 << "\n\n"
            << "h2 = " << h2 << '\n'
            << "Rectangle Method: " << R2 << '\n'
            << "Trapezoid Method: " << T2 << '\n'
            << "Simpson Method: " << S2 << "\n\n"
            << "Runge Romberg Richardson Method:\n" 
            << "for Rectangle Method: " << posteriori_R << '\n'
            << "for Trapezoid Method: " << posteriori_T << '\n'
            << "for Simpson Method: " << posteriori_S << '\n';
}

int main() {
    task1();
    task2();
    task3();
    task4();
    task5();
}