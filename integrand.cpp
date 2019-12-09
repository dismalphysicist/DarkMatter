#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp> 
#include <valarray>
#include "integrand.h"

using namespace std;

Integrand::Integrand (double m_, double T_, double V_, double A_, double Vtil_, double Atil_) {
    m = m_;
    T = T_;
    v0 = k_B*T/m;
    V = V_;
    A = A_;
    Vtil = Vtil_;
    Atil = Atil_;
}

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
    double factor0 = 1/(sqrt(2*pi) * pow(m,4) * pow(v0,3)); //lois flower 
    //double factor0 = 1.0;

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
    // Nonrelativistic 
    double velocity = sqrt(s - 4.0*(m*m))/m; 
    double prefactor = sqrt(2/pi)/(2*pow(v0,3)*pow(m,3)); 
	double weight = sqrt(s - 4.0*(m*m)) * exp(-(s-4.0*m*m)/(2.0*m*m*v0*v0));

    // Relativistic 
    // double velocity = sqrt(s*(s-4*m*m))/(s-2*m*m); 
    // double prefactor = 1/(8*T*pow(m,3)*pow(boost::math::cyl_bessel_k(2,m/T),2));
    // double weight = sqrt(s - 4.0*(m*m))*(s-2*m*m)*boost::math::cyl_bessel_k(1,sqrt(s)/T); 

    return prefactor*sigma(s)*velocity*weight;
}

valarray<double> Integrand::sigmav (valarray<double> s) {
    valarray<double> ret = valarray<double>(s.size());
    
    for (int i = 0; i < s.size(); i++) {
        ret[i] = Integrand::sigmav(s[i]);
    }
    return ret;
}

/*
int main() {
    Integrand inte(1,1,1,1,1,1); 

    valarray<double> test = {4,5,6,7,8};
    for(double i : test) cout << i << ", ";
    cout << endl;

    for(double i : test) cout << inte.sigmav(i) << ", "; //element by element
    cout << endl;

    for(double a : inte.sigmav(test)) cout << a << ", "; //whole array
    cout << endl;
    return 0;
}
*/