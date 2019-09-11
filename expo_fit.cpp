#include <iostream>
#include <cmath>
#include <valarray>
#include "expo_fit.h"

using namespace std;

double Expo_fit::g1 (double s) {
    if (s < peak_0) {  
        return height * exp(k1 * (s - peak_0));
    }
    else {
        return height * exp(k2 * (peak_0 - s));
    }
}

valarray<double> Expo_fit::g1 (valarray<double> ss) {
    /* g(s) implemented here as two connected exponential functions with the same peak at s=peak_0. */
    valarray<double> ret = valarray<double>(ss.size());
    
    for (int i = 0; i < ss.size(); i++) {
        ret[i] = Expo_fit::g1(ss[i]); 
    }
    return ret;
}

double Expo_fit::rho1 (double s) {
    if (s < peak_0) {  
        return height/k1 * exp(k1*(s-peak_0));
    }
    else {
        return height/k1 - height/k2 * exp(k2*(peak_0-s)) + height/k2;
    }
}

valarray<double> Expo_fit::rho1 (valarray<double> ss) {
    /* rho(s). It is essentially indefinite integral of g(s) */
    valarray<double> ret = valarray<double>(ss.size());
    
    for (int i = 0; i < ss.size(); i++) {
        ret[i] = Expo_fit::rho1(ss[i]); 
    }
    return ret;
}

double Expo_fit::s1 (double rho) {
    if (rho < height/k1) {  
        return log(k1*rho/height)/k1 + peak_0;
    }
    else {
        return -log(1 + k2/k1 - k2*rho/height)/k2 + peak_0;
    }
}

valarray<double> Expo_fit::s1 (valarray<double> rhos) {
    /* s(rho). The inverse function of rho(s). It is used to generate importance sampling on s. */
    valarray<double> ret = valarray<double>(rhos.size());
    
    for (int i = 0; i < rhos.size(); i++) {
        ret[i] = Expo_fit::s1(rhos[i]); 
    }
    return ret;
}

double Expo_fit::g2 (double s) {
    return height_2 / (b1 * (s-peak_2)*(s-peak_2) + 1.0);
}
valarray<double> Expo_fit::g2 (valarray<double> ss) {
    valarray<double> ret = valarray<double>(ss.size());
    
    for (int i = 0; i < ss.size(); i++) {
        ret[i] = Expo_fit::g2(ss[i]); 
    }
    return ret;
}

double Expo_fit::rho2 (double s) {
    return height_2/sqrt(b1) * atan((s-peak_2)*sqrt(b1));
}

valarray<double> Expo_fit::rho2 (valarray<double> ss) {
    valarray<double> ret = valarray<double>(ss.size());
    
    for (int i = 0; i < ss.size(); i++) {
        ret[i] = Expo_fit::rho2(ss[i]); 
    }
    return ret;
}
double Expo_fit::s2 (double rho) {
    return peak_2 + 1.0/sqrt(b1) * tan(sqrt(b1) * rho / height_2);
}

valarray<double> Expo_fit::s2 (valarray<double> rhos) {
    valarray<double> ret = valarray<double>(rhos.size());
    
    for (int i = 0; i < rhos.size(); i++) {
        ret[i] = Expo_fit::s2(rhos[i]); 
    }
    return ret;
}

double Expo_fit::g_tot (double s, double alpha1, double alpha2, double norm_factor_1, double norm_factor_2) {
    return alpha1*norm_factor_1*g1(s) + alpha2*norm_factor_2*g2(s);
}

valarray<double> Expo_fit::g_tot (valarray<double> ss, double alpha1, double alpha2, double norm_factor_1, double norm_factor_2) {
    valarray<double> ret = valarray<double>(ss.size());
    
    for (int i = 0; i < ss.size(); i++) {
        ret[i] = Expo_fit::g_tot(ss[i], alpha1, alpha2, norm_factor_1, norm_factor_2); 
    }
    return ret;
}

/*
int main() {
    Expo_fit e;
    double test_s = {20};
    double rho = e.rho2(test_s);
    cout << rho << endl;
    
    double s = e.s2(rho);
    cout << s << endl;
    cout << e.g_tot(5,0.5,0.5,2,2) << endl; 
    return 0;
}
*/