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
    double peak_2 = 30.0;
    double height_2 = 0.00002;
    double b1 = 7.0;

    public: 

    double g1 (double s);
    valarray<double> g1 (valarray<double> ss);

    double rho1 (double s);
    valarray<double> rho1 (valarray<double> ss);

    double s1 (double rho);
    valarray<double> s1 (valarray<double> rhos);


    double g2 (double s);
    valarray<double> g2 (valarray<double> ss);

    double rho2 (double s);
    valarray<double> rho2 (valarray<double> ss);

    double s2 (double rho);
    valarray<double> s2 (valarray<double> rhos);

    double g_tot (double s, double alpha1, double alpha2, double norm_factor_1, double norm_factor_2);
    valarray<double> g_tot (valarray<double> ss, double alpha1, double alpha2, double norm_factor_1, double norm_factor_2);
}; 