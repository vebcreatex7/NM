#include "calculation.hpp"




void Print1(std::vector<ld> const & x, std::vector<ld> const y, std::vector<ld> const & est, std::ostream & myFile) {
    for (size_t i = 0; i != x.size(); i++)
        myFile <<  std::setw(7)  << i << ' ';
    myFile << '\n';

    myFile << std::setw(6) << std::left << "x_k:";
    for (size_t i = 0; i != x.size(); i++)
        myFile << std::setprecision(5) << std::fixed << x[i] << ' ';
    myFile << '\n';
    myFile << std::setw(6) << std::left << "y_k:";
    for (size_t i = 0; i != x.size(); i++)
        myFile  << y[i] << ' ';
    myFile << '\n';

    myFile << std::setw(6) << std::left << "eps_k:";
    for (size_t i = 0; i != x.size(); i++)
        myFile  << eps(y[i], y_exact1(x[i])) << ' ';
    myFile << '\n';
    

    myFile << std::setw(6) << std::left << "Runge:";
    for (size_t i = 0; i != x.size(); i++)
        myFile  << est[i] << ' ';
    myFile << '\n';
}


void Task1() {
    std::fstream myFile("../tasks/t1.txt", std::ios::in);
    ld x_0, x_1, y_0, z_0;
    myFile >> x_0 >> x_1 >> y_0 >> z_0;
    ld h;
    myFile >> h;
    myFile.close();



    //Euler
    std::ofstream log("../logs/log1.1.txt");
    auto [x1, y1] = Euler_Method(h, x_0, x_1, y_0, z_0, &log);
    auto Runge_Euler = Runge_Estimate_for_Euler(y1, h, x_0, x_1, y_0, z_0);
    log.close();
    myFile.open("../ans/a1.txt", std::ios::out);
    myFile << "Euler's Method\n";
    Print1(x1, y1, Runge_Euler, myFile);
    

    //Runge-Kutta
    log.open("../logs/log1.2.txt");
    auto [x2, y2, other] = Runge_Kutta_Method(h, x_0, x_1, y_0, z_0, &log);
    auto Runge_Runge_Kutta = Runge_Estimate_for_Runge_Kutta(y2, h, x_0, x_1, y_0, z_0);
    log.close();
    myFile << "\n\nRunge-Kutta's Method\n";
    Print1(x2, y2, Runge_Runge_Kutta, myFile);
    

    //Adams
    log.open("../logs/log1.3.txt");
    auto [x3, y3] = Adams_Method(h, x_0, x_1, y_0, z_0, &log);
    auto Runge_Adams = Runge_Estimate_for_Adams(y3, h, x_0, x_1, y_0, z_0);
    log.close();
    myFile << "\n\nAdams's Method\n";
    Print1(x3, y3, Runge_Adams, myFile);

    myFile.close();
    
}



//Task2

void Print2(std::vector<ld> const & x, std::vector<ld> const y, std::vector<ld> const & est, std::ostream & myFile) {
    for (size_t i = 0; i != x.size(); i++)
        myFile <<  std::setw(7)  << i << ' ';
    myFile << '\n';

    myFile << std::setw(6) << std::left << "x_k:";
    for (size_t i = 0; i != x.size(); i++)
        myFile << std::setprecision(5) << std::fixed << x[i] << ' ';
    myFile << '\n';
    myFile << std::setw(6) << std::left << "y_k:";
    for (size_t i = 0; i != x.size(); i++)
        myFile  << y[i] << ' ';
    myFile << '\n';

    myFile << std::setw(6) << std::left << "eps_k:";
    for (size_t i = 0; i != x.size(); i++)
        myFile  << eps(y[i], y_exact2(x[i])) << ' ';
    myFile << '\n';
    

    myFile << std::setw(6) << std::left << "Runge:";
    for (size_t i = 0; i != x.size(); i++)
        myFile  << est[i] << ' ';
    myFile << '\n';
}

void Task2() {
    std::fstream myFile("../tasks/t2.txt", std::ios::in);
    ld x_0, x_1;
    ld a_0, a_1, alpha;
    ld b_0, b_1, beta;
    myFile >> x_0 >> x_1;
    myFile >> a_0 >> a_1 >> alpha;
    myFile >> b_0 >> b_1 >> beta;
    ld h;
    myFile >> h;
    myFile.close();

    std::ofstream log("../logs/log2.1.txt");
    auto [x1, y1] = Finite_Difference_Method(h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta, &log);
    auto Runge1 = Runge_Estimate_for_Finite_Defference(y1, h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta);
    log.close();
    myFile.open("../ans/a2.txt", std::ios::out);
    myFile << "Finite Defference Method\n";
    Print2(x1, y1, Runge1, myFile);

    log.open("../logs/log2.2.txt");
    auto [x2, y2] = Shooting_Method(h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta, &log);
    auto Runge2 = Runge_Estimate_for_Shooting(y2, h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta);
    log.close();
    myFile << "\n\nShooting Method\n";
    Print2(x2, y2, Runge2, myFile);

    myFile.close();
}


//Shooting






int main() {


    Task1();
    Task2();
    //Shooting();

    return 0;
}