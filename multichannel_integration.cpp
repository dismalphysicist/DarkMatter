#include <iostream>
#include <cmath>
#include <valarray>
#include <random>
#include <ctime>
#include "integrand.h"
#include "expo_fit.h"

using namespace std;

const double pi = acos(-1);

//////////// DARK MATTER PARAMETERS /////////////
double m = 20; //dark matter mass 
double T = 2; //temperature
double V = -0.1; //DM vector coupling constant 0.1 * 1
double Vtil = sqrt(4*pi/128) * (-0.5 + 2*0.23)/ (2*sqrt(0.23)*0.88); //electron vector coupling constant 
double A = 0.1; //DM axial c.c. 0.1 * 1
double Atil = sqrt(4*pi/128) *0.5 / (2*0.23*0.88); //electron axial c.c. 

//////////// INTEGRATION PARAMETERS /////////////
double height = 1.5e-11; //height of exponentials modelling velocity peak 
double k2 = 0.25e-3; //exp(k1 s)
double k1 = 0.5e-3; //exp(-k2 s)

double height_2 = 4e-6; //height of g2(s) modelling cross section peak 

double start_point = 4*m*m + 1; //cutoff minimum energy, so sigmav does not blow up or become imaginary
double end_point = 10000; //a good approximate for infinity, since the tail is very small beyond this point 

double alpha_1 = 0.5; //probability of sampling from g1(s)
double alpha_2 = 0.5; //probability of sampling from g2(s)

//these are initialised by the initialise() function 
double int1_1, int2_1, area_1, norm_factor_1; 
double int1_2, int2_2, area_2, norm_factor_2;

void initialise(Integrand inte, Expo_fit fit,
double &int1_1, double &int2_1, double &int1_2, double &int2_2,
double &area_1, double &area_2, double &norm_factor_1, double &norm_factor_2,
double start = start_point, double end = end_point) {
    /* Initialisation: specify range of integration. */ 

    int1_1 = fit.rho1(start);
    int2_1 = fit.rho1(end);
    area_1 = int2_1 - int1_1; 
    norm_factor_1 = 1.0/area_1;

    int1_2 = fit.rho2(start);
    int2_2 = fit.rho2(end);
    area_2 = int2_2 - int1_2; 
    norm_factor_2 = 1.0/area_2;
}

double monteCarloIntegrate(int N, int rounds, double &err, double start=start_point, double end=end_point) {
    /* Monte Carlo integrator by importance sampling
    N: number of points per round
    rounds: number of rounds
    start, end: floats, range of integration */ 

    //passing relevant parameters 
    Integrand inte(m,T,V,A,Vtil,Atil);
    Expo_fit fit(m,T,height,k1,k2,inte.M,height_2,inte.gamma);

    initialise(inte, fit, int1_1,int2_1,int1_2,int2_2,area_1,area_2,norm_factor_1,norm_factor_2);
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
        valarray<double> rhos2 = valarray<double>(N);
        for (int i=0; i<N; i++) {
            rhos1[i] = uni1(rng);
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
            ys[i] = inte.sigmav(ss[i])/fit.g_tot(ss[i], alpha_1, alpha_2, norm_factor_1, norm_factor_2);
        }
        /*
        for (int i=0;i<N;i++) {
            cout << ys[i] << ", ";
        }
        cout << endl;
        */

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
    double mc = monteCarloIntegrate(1000000,100,err);

    cout << mc << " GeV^-2" << endl;
    cout << err << " GeV^-2" << endl;

    cout << mc*0.389/1e-9 << " picobarns" << endl;
    cout << err*0.389/1e-9 << " picobarns" << endl;
    return 0;
}