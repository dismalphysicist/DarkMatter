#include <valarray>
#include <cmath>

using namespace std;

class Integrand {
    const double pi = acos(-1);

    double m; //dark matter mass 
    double T; //temperature
    double k_B = 1; //boltzmann constant 
    double v0;
    double V;
    double A;
    double Vtil;
    double Atil;


    public: 
    double M = 91.2; //Z boson mass 
    double gamma = 2.50; //Z boson width 
    
    Integrand (double m_, double T_, double V_, double A_, double Vtil_, double Atil_);

    double f (double s);
    valarray<double> f (valarray<double> s);
    double f_multi (double s);
    valarray<double> f_multi (valarray<double> s);

    double sigma (double s);
    valarray<double> sigma (valarray<double> s);

    double sigmav (double s);
    valarray<double> sigmav (valarray<double> s);
};
