#include <valarray>

using namespace std;

class Integrand {
    const double pi = acos(-1);

    public: 
    //double f (double s);
    valarray<double> f (valarray<double> s);
    //double f_multi (double s);
    valarray<double> f_multi (valarray<double> s);
};