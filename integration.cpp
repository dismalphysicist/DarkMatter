#include <iostream>
#include <cmath>
#include <valarray>
#include <random>
#include <ctime>
#include "integrand.h"
#include "expo_fit.h"

using namespace std;

const double pi = acos(-1);

double m = 20; //dark matter mass 
double T = 20; //temperature
double V = 0.1 * 1;
double Vtil = sqrt(4*pi/128) * (1 + 2*0.23)/ (2*sqrt(0.23)*0.88);
double A = 0.1 * 1;
double Atil = sqrt(4*pi/128) / (2*sqrt(0.23)*0.88);

double start_point = 4.0;
double end_point = 25.0;

double monteCarloIntegrate(int N, int rounds, double &err, double start=start_point, double end=end_point) {
    double ans = 0;
    double sum_y = 0;
    double sum_ysq = 0;

    Integrand inte(m,T,V,A,Vtil,Atil);

    Expo_fit fit(m,1,1,1,1,0,0,0);

    double area = fit.rho1(end) - fit.rho1(start);
    double norm_factor = 1.0/area;
    cout << "Normalisation factor = " << norm_factor << endl;

    //setting up random number generator 
    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    uniform_real_distribution<double> uni(fit.rho1(start),fit.rho1(end)); // guaranteed unbiased

    for (int i = 0; i < rounds; i++) {
        
        //initialising array of rhos with random numbers 
        valarray<double> rhos = valarray<double>(N);
        for (int i=0; i<N; i++) {
            rhos[i] = uni(rng);
        }

        //getting suitably distributed s values from uniform rhos 
        valarray<double> ss = fit.s1(rhos);

        //y = f/g with normalisation factor 
        valarray<double> ys = inte.f(ss)/(norm_factor * fit.g1(ss));

        //summing up arrays, keeping running total
        sum_y += ys.sum();
        sum_ysq += (ys*ys).sum();

        //result and error 
        ans = sum_y / ((i+1)*N);
        err = sqrt((sum_ysq/((i+1)*N) - (sum_y/((i+1)*N))*(sum_y/((i+1)*N)))/((i+1)*N));
    }
    return ans;
}

int main() {
    //cout << f(5) << endl;
    double err;

    clock_t begin = clock();

    cout << monteCarloIntegrate(1000000,100,err) << endl;

    clock_t end = clock();

    cout << err << endl;
    cout << double(end-begin)/CLOCKS_PER_SEC << " s" << endl;
    return 0;
}