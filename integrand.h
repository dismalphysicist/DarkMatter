#include <valarray>
#include <cmath>

using namespace std;

class Integrand {
    const double pi = acos(-1);
    double m = 20; //dark matter mass 
    double T = 20; //temperature
    double k_B = 1; //boltzmann constant 
    double v0 = k_B*T/m;
    double M = 91.2; //Z boson mass 
    double gamma = 2.50; //Z boson width 
    double V = 1;
    double Vtil = 1;
    double A = 0.5;
    double Atil = 0.5;


    public: 
    double f (double s);
    valarray<double> f (valarray<double> s);
    double f_multi (double s);
    valarray<double> f_multi (valarray<double> s);

    double sigma (double s);
    valarray<double> sigma (valarray<double> s);

    double sigmav (double s);
    valarray<double> sigmav (valarray<double> s);
};