#include <iostream>
#include <cmath>
#include <valarray>

using namespace std;

class Expo_fit {
    private:
    double peak_0;
    double height;
    double k1 = 0.5e-3;
    double k2 = 0.25e-3;

    double peak_2;
    double height_2;
    double width_term;

    public: 
    Expo_fit (double m, double height_, double k1_, double k2_, double M, double height_2_, double width);

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