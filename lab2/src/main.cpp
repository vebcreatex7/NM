#include "calculation.hpp"
#include "matrix.hpp"




void Task1() {

    long double eps;
    std::fstream myFile("../tests/t1.txt", std::ios::in);


    myFile >> eps;
    myFile.close();


    std::fstream log1("../logs/log1.1.txt", std::ios::out);
    std::fstream log2("../logs/log1.2.txt", std::ios::out);
    auto [x1, count1] = Newton_Method(0, 0.5, eps, log1);
    auto [x2, count2] = Simple_Iterations_Method(0, 1, eps, log2);
    log1.close();
    log2.close();

    myFile.open("../ans/a1.txt");
    myFile << "Newton Method:\n" << "x = " << x1 << "\nNumber of iterations: " << count1 << "\n\n";
    myFile << "Simple Iterations Method:\n" << "x = " << x2 << "\nNumber of iterations: " << count2 << "\n";
    
    myFile.close();
}

void Task2() {
    long double eps;
    std::fstream myFile("../tests/t2.txt", std::ios::in);

    myFile >> eps;
    myFile.close();

    std::fstream log1("../logs/log2.1.txt", std::ios::out);
    std::fstream log2("../logs/log2.2.txt", std::ios::out);
    auto [x1, count1] = Dimentional_Newton_Method(eps, log1);
    auto [x2, count2] = Dimentional_Simple_Iterations_Method(eps, log2);
    log1.close();
    log2.close();

    myFile.open("../ans/a2.txt");
    myFile << "Newton Method:\n" << "x:\n" << x1 << "\nNumber of iterations: " << count1 << "\n\n"
            << "Simple Iterations Method:\n" << "x:\n" << x2 << "\nNumber of iterations: " << count2 << "\n\n";

}
int main() {
    Task1();
    Task2();

    return 0;
}