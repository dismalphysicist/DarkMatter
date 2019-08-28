#include <iostream>
#include <cmath>
#include <valarray>

using namespace std;

class Expo_fit {
    private:
    double peak0 = 2 + 2*sqrt(3);
    double height_ = 0.1295;
    double k2_ = 0.3;
    double k1_ = 1;

    public: 
    valarray<double> g1 (valarray<double> s) {
        /* g(s) implemented here as two connected exponential functions with the same peak at peak0. */
        valarray<double> ret = valarray<double>(s.size());
        
        for(double i : ret) cout << i << ", "; //debugging 
        cout << endl; //debugging 
    }
}