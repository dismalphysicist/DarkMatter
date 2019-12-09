#include <valarray>
#include <cmath>
#ifndef INTEGRAND
#define INTEGRAND 
#include "integrand.h"
#endif

using namespace std;

class IntegrandY {
    private: 

    const double pi = acos(-1);

    double m; //dark matter mass 
    double T; //temperature
    double k_B = 1; //boltzmann constant 
    double v0;

    public: 

    IntegrandY (double m_, double T_);

    double boost (double y);
    valarray<double> boost (valarray<double> y);

    inline double sqr(double x) {return x*x;}
};