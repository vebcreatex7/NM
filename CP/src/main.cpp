#include "calculation.hpp"

void Print(std::vector<ld> const& x, std::vector<ld> const & y, std::ostream & myFile) {
    for (auto a : x) 
        myFile << a << ' ';
    myFile << '\n';

    for (auto a : y)
        myFile << a << ' ';
    myFile << '\n';
    
}

void solve() {
    std::fstream myFile("../input.txt", std::ios::in);
    ld x_0, x_1;
    ld a_0, a_1, alpha;
    ld b_0, b_1, beta;
    myFile >> x_0 >> x_1;
    myFile >> a_0 >> a_1 >> alpha;
    myFile >> b_0 >> b_1 >> beta;
    ld h;
    myFile >> h;
    myFile.close();

    std::ofstream log("../logs/log1.txt");
    auto [x1, y1] = Shooting_Method(h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta, &log);
    log.close();

     log.open("../logs/log2.txt");
    auto [x2, y2] = Finite_Difference_Method(h, x_0, x_1, a_0, a_1, alpha, b_0, b_1, beta, &log);
    log.close();

    myFile.open("../ans.txt", std::ios::out);
    myFile << "Finite Difference Method:\n";
    Print(x2, y2, myFile);
    myFile << "\n\nShooting Method:\n";
    Print(x1, y1, myFile);
    myFile.close();

}


int main() {
    solve();
}