#include <cmath>
#include <fstream>

class TNewton_Method {

public:
    TNewton_Method(long double, long double);
    ~TNewton_Method();
    long double f(long double x);
    long double First_Derivative(long double x);
    long double Second_Derivative(long double x);
    long double Solution();
    
private:
    long double a, b;

};