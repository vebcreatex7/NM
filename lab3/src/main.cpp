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

int main() {
    //task1();
    task2();
}