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

/*
int main() {
    cout << log(1) << endl;
    cout << log(2.71828) << endl;
    Expo_fit e;
    valarray<double> test = {0,5,10,15};
    e.g1(test);
    
    valarray<double> test_rhos = {0.1,0.2,0.3,0.4,0.5};
    valarray<double> ss = e.s1(test_rhos);
    for(double i : ss) cout << i << ", ";
    cout << endl; 
    return 0;
}
*/