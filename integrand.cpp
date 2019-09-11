#include <iostream>
#include <cmath>
#include <valarray>
#include "integrand.h"

using namespace std;


double Integrand::f (double s) {

    //constant scalar coefficient
    //double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    double factor0 = 1.0;

    double factor1 = (s - 4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    double factor2 = 1.0/s;

    return factor0 * factor1 * factor2;   
}


valarray<double> Integrand::f (valarray<double> s) {

    //constant scalar coefficient
    //double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4)); //liu zhenning 
    //double factor0 = 1/(sqrt(2*pi) * pow(m,4) * pow(v0,3)); //lois flower 
    double factor0 = 1.0;

    valarray<double> factor1 = (s-4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    valarray<double> factor2 = 1.0/s;

    return factor0 * factor1 * factor2;   
}


double Integrand::f_multi (double s) {

    //constant scalar coefficient
    double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    //double factor0 = 1.0;

    double factor1 = (s - 4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    double factor2 = 1.0 / (s * ((s-30.0)*(s-30.0) + 0.1));

    return factor0 * factor1 * factor2;   
}


valarray<double> Integrand::f_multi (valarray<double> s) {

    //constant scalar coefficient
    double factor0 = sqrt(2/pi) / (pow(v0,3)*pow(m,4));
    //double factor0 = 1.0;

    valarray<double> factor1 = (s-4*m*m) * exp(-(s-4*m*m)/(2*m*m * v0*v0)); //velocity and distribution

    valarray<double> factor2 = 1.0 / (s * ((s-30.0)*(s-30.0) + 0.1));

    return factor0 * factor1 * factor2;   
}

double Integrand::sigma (double s) {
    double prefactor = 1/(8*pi);
    double factor1 = 1/sqrt(1-4*m*m/s);
    double factor2 = 1/((s-M*M)*(s-M*M) + M*M * gamma*gamma);
    double part1 = (V*V+A*A)*(Vtil*Vtil+Atil*Atil)*(s/3 - m*m/3);
    double part2 = (V*V-A*A)*(Vtil*Vtil+Atil*Atil)*2*m*m; 
    
    return prefactor*factor1*factor2*(part1+part2);
}

valarray<double> Integrand::sigma (valarray<double> s) {
    valarray<double> ret = valarray<double>(s.size());
    
    for (int i = 0; i < s.size(); i++) {
        ret[i] = Integrand::sigma(s[i]);
    }
    return ret;
}

double Integrand::sigmav (double s) {
    double prefactor = 1/(sqrt(2*pi) * pow(m,4) * pow(v0,3));
    double velocityfactor = (s - 4.0*(m*m)) * exp(-(s-4.0*m*m)/(2.0*m*m*v0*v0));
    return prefactor*sigma(s)*velocityfactor;
}

valarray<double> Integrand::sigmav (valarray<double> s) {
    double prefactor = 1/(sqrt(2*pi) * pow(m,4) * pow(v0,3));
    valarray<double> velocityfactor = (s - 4.0*(m*m)) * exp(-(s-4.0*m*m)/(2.0*m*m*v0*v0));
    return prefactor*sigma(s)*velocityfactor;
}

/*
int main() {
    Integrand inte; 

    valarray<double> test = {4,5,6,7,8};
    for(double i : test) cout << i << ", ";
    cout << endl;

    for(double i : test) cout << inte.f_multi(i) << ", "; //element by element
    cout << endl;

    for(double a : inte.f_multi(test)) cout << a << ", "; //whole array
    cout << endl;
    return 0;
}
*/