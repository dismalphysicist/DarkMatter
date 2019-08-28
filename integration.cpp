#include <iostream>
#include <cmath>
#include <valarray>
#include <random>
#include <time.h>
#include "integrand.h"
#include "expo_fit.h"

using namespace std;

double start_point = 4.0;
double end_point = 25.0;

double monteCarloIntegrate(int N, int rounds, double &err, double start=start_point, double end=end_point) {
    double ans = 0;
    double sum_y = 0;
    double sum_ysq = 0;
    double s,y;

    srand(time(NULL));

    Expo_fit fit;
    double area = fit.rho1(end) - fit.rho1(start);
    double norm_factor = 1.0/area;
    cout << "Normalisation factor = " << norm_factor << endl;

    for (int i = 0; i < rounds; i++) {
        //initialising array of rhos with random numbers 
        valarray<double> rhos = valarray<double>(N);

        random_device rd;     // only used once to initialise (seed) engine
        mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
        uniform_real_distribution<double> uni(fit.rho1(start),fit.rho1(end)); // guaranteed unbiased

        for (int i=0; i<N; i++) {
            rhos[i] = uni(rng);
        }
        for(double i : rhos) cout << i << ", ";
        cout << endl;

        //getting suitably distributed s values from uniform rhos 
        valarray<double> ss = fit.s1(rhos);
        for(double i : ss) cout << i << ", ";
        cout << endl;

        //y = f/g with normalisation factor 
        valarray<double> ys = f(ss)/(norm_factor * fit.g1(ss));
        for(double i : ys) cout << i << ", ";
        cout << endl;

        //summing up arrays, keeping running total
        sum_y += ys.sum();
        valarray<double> ysq = ys*ys;
        sum_ysq += ysq.sum();

        //result and error 
        ans = sum_y / ((i+1)*N);
        err = sqrt(sum_ysq/((i+1)*N) - pow(sum_y/((i+1)*N),2))/((i+1)*N);
    }
    return ans;
}

int main() {
    cout << f(5) << endl;
    double err;
    cout << monteCarloIntegrate(100,1,err) << endl;
    cout << err << endl;
    return 0;
}