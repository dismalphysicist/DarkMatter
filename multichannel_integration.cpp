#include <iostream>
#include <cmath>
#include <valarray>
#include <random>
#include <ctime>
#include "integrand.h"
#include "expo_fit.h"

using namespace std;

double start_point = 4.0;
double end_point = 40.0;

double alpha_1 = 0.5;
double alpha_2 = 0.5;

double int1_1=1.0, int2_1=1.0, area_1=1.0, norm_factor_1=1.0; 
double int1_2=1.0, int2_2=1.0, area_2=1.0, norm_factor_2=1.0;

void initialise(Integrand inte, Expo_fit fit, double start = start_point, double end = end_point) {
    /* Initialisation: specify range of integration. */ 

    int1_1 = fit.rho1(start);
    int2_1 = fit.rho1(end);
    double area_1 = int2_1 - int1_1; 
    double norm_factor_1 = 1.0/area_1;

    int1_2 = fit.rho2(start);
    int2_2 = fit.rho2(end);
    double area_2 = int2_2 - int1_2; 
    double norm_factor_2 = 1.0/area_2;
}

double monteCarloIntegrate(int N, int rounds, double &err, double start=start_point, double end=end_point) {
    /* Monte Carlo integrator by importance sampling
    N: number of points per round
    rounds: number of rounds
    fit_1: the first fitting function
    fit_2: the second fitting function
    inte: function, the integrand
    start, end: floats, range of integration */ 

    Integrand inte;
    Expo_fit fit;

    initialise(inte, fit);
    //cout << "Initialised" << endl; //debugging 

    double ans = 0;
    double sum_y = 0;
    double sum_ysq = 0;
    double s,y;

    //setting up random number generator 
    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    uniform_real_distribution<double> uni1 (int1_1,int2_1); // guaranteed unbiased
    uniform_real_distribution<double> uni2 (int1_2,int2_2); 
    uniform_real_distribution<double> rands (0,1); 
    //cout << "Random generators set up" << endl; //debugging 

    for (int j = 0; j < rounds; j++) {

        //initialising channels 
        valarray<double> channels = valarray<double>(N);
        for (int i=0; i<N; i++) {
            channels[i] = rands(rng);
        }

        //initialising arrays of rhos with random numbers 
        valarray<double> rhos1 = valarray<double>(N);
        for (int i=0; i<N; i++) {
            rhos1[i] = uni1(rng);
        }

        valarray<double> rhos2 = valarray<double>(N);
        for (int i=0; i<N; i++) {
            rhos2[i] = uni2(rng);
        }

        //getting suitably distributed s values from uniform rhos 
        valarray<double> ss = valarray<double>(N);
        valarray<double> ys = valarray<double>(N);

        for (int i=0; i<N; i++) {
            if (channels[i] < alpha_1) {
                ss[i] = fit.s1(rhos1[i]);
            }
            else { 
                ss[i] = fit.s2(rhos2[i]);
            }

            //y = f/g with normalisation factor
            ys[i] = inte.f_multi(ss[i])/fit.g_tot(ss[i], alpha_1, alpha_2, norm_factor_1, norm_factor_2);
        }

        sum_y += ys.sum();
        sum_ysq += (ys*ys).sum();

        //result and error 
        ans = sum_y / ((j+1)*N);
        err = sqrt((sum_ysq/((j+1)*N) - (sum_y/((j+1)*N))*(sum_y/((j+1)*N)))/((j+1)*N));
    
        //weights? not sure what this all does
        /*
        double dw1 = (fit.g1(ss) *ys*ys * norm_factor_1 / fit.g_tot(ss[i], alpha_1, alpha_2, norm_factor_1, norm_factor_2)).sum();
        double dw2 = (fit.g2(ss) *ys*ys * norm_factor_1 / fit.g_tot(ss[i], alpha_1, alpha_2, norm_factor_1, norm_factor_2)).sum();
    
        double alpha1n = alpha_1 * sqrt(dw1);
        double alpha2n = alpha_2 * sqrt(dw2);

        double alpha1 = alpha1n / (alpha1n+alpha2n);                         
        double alpha2 = alpha2n / (alpha1n+alpha2n);     
        */
    }
    return ans;
}

int main() {
    double err;
    cout << monteCarloIntegrate(100,1,err) << endl;
    cout << err << endl;
    return 0;
}