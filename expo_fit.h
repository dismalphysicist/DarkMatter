#include <iostream>
#include <cmath>
#include <valarray>

using namespace std;

class Expo_fit {
    private:
    double peak_0 = 2 + 2*sqrt(3.0);
    double height = 0.1295;
    double k2 = 0.3;
    double k1 = 1;

    public: 
    double g1 (double s);
    valarray<double> g1 (valarray<double> ss);
    double rho1 (double s);
    valarray<double> rho1 (valarray<double> ss);
    double s1 (double rho);
    valarray<double> s1 (valarray<double> rhos);
}; 