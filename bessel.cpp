#include <iostream> 
#include <cmath>
#include <valarray>

using namespace std;

double pi = acos(-1);

double cyl_bessel_k (int v, double z) {
	return sqrt(pi/(2*z)) * exp(-z); //large z approximation - depends on v 
}

int main() {
	valarray<double> zs = {2,4,6,8,10,12,14};
	valarray<double> bessels1 = valarray<double>(zs.size());

	for (int i=0; i<zs.size(); i++) {
		bessels1[i] = cyl_bessel_k(1,zs[i]); 
		cout << bessels1[i] << ", ";
	}
	cout << endl; 

	return 0;
}
