#include <iostream>
#include <cmath>
#include <valarray>
#include "integrand_y.h" 

using namespace std;

IntegrandY::IntegrandY (double m_, double T_) {
    m = m_;
    T = T_;
    v0 = k_B*T/m;
}

double IntegrandY::boost (double y) {
    return exp(-2/sqr(v0) * sqr(tanh(y))) * sqr(tanh(y)) * (1+sqr(tanh(y))); 
}

valarray<double> IntegrandY::boost (valarray<double> y) {
    valarray<double> ret = valarray<double>(y.size());
    
    for (int i = 0; i < y.size(); i++) {
        ret[i] = IntegrandY::boost(y[i]);
    }
    return ret;
}

// int main() {
//     IntegrandY inte = IntegrandY(10,2); 

//     valarray<double> test = {-1e-4,-1e-2,0,1e-2,1e-4};

//     for(double i : test) cout << inte.boost(i) << ", "; //element by element
//     cout << endl;

//     for(double a : inte.boost(test)) cout << a << ", "; //whole array
//     cout << endl;
//     return 0;
// }