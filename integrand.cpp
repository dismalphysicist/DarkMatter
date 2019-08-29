#include <iostream>
#include <cmath>
#include <valarray>
#include "integrand.h"

using namespace std;

/*
double Integrand::f (double s) {

    double m = 1;
    double v0 = 1;

    //constant scalar coefficient
    //double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    double factor0 = 1.0;

    double factor1 = (s - 4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    double factor2 = 1.0/s;

    return factor0 * factor1 * factor2;   
}
*/

valarray<double> Integrand::f (valarray<double> s) {

    double m = 1;
    double v0 = 1;

    //constant scalar coefficient
    //double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    double factor0 = 1.0;

    valarray<double> factor1 = (s-4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    valarray<double> factor2 = 1.0/s;

    return factor0 * factor1 * factor2;   
}

/*
double Integrand::f_multi (double s) {
    double m = 1;
    double v0 = 1;

    //constant scalar coefficient
    //double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    double factor0 = 1.0;

    double factor1 = (s - 4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    double factor2 = 1.0 / (s * ((s-30.0)*(s-30.0) + 0.1));

    return factor0 * factor1 * factor2;   
}
*/

valarray<double> Integrand::f_multi (valarray<double> s) {
    double m = 1;
    double v0 = 1;

    //constant scalar coefficient
    //double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    double factor0 = 1.0;

    valarray<double> factor1 = (s-4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    valarray<double> factor2 = 1.0 / (s * ((s-30.0)*(s-30.0) + 0.1));

    return factor0 * factor1 * factor2;   
}

/*
int main() {
    valarray<double> test = {4,5,6,7,8};
    for(double i : test) cout << i << ", ";
    cout << endl;

    for(double i : test) cout << f(i) << ", "; //element by element
    cout << endl;

    for(double a : f(test)) cout << a << ", "; //whole array
    cout << endl;
    return 0;
}
*/