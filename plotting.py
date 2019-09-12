import numpy as np 
from matplotlib import pyplot 

#parameters
m = 20.0 #dark matter mass 
T = 20.0 #temperature
k_B = 1 #boltzmann constant 
v0 = np.sqrt(k_B*T/m) 
M = 91.2 #Z boson mass 
gamma = 2.50 #Z boson width 
V = 1.0
Vtil = 1.0 
A = 0.5
Atil = 0.5


#expression for cross section 
def sigma (s):
    prefactor = 1/(8*np.pi)
    factor1 = 1/np.sqrt(1-4*m**2/s)
    factor2 = 1/((s-M**2)**2 + M**2 * gamma**2)
    part1 = (V**2+A**2)*(Vtil**2+Atil**2)*(s/3 - m**2/3)
    part2 = (V**2-A**2)*(Vtil**2+Atil**2)*2*m**2 
    
    return prefactor*factor1*factor2*(part1+part2)

#cross section x velocity distribution 
def sigmav (s):
    prefactor = 1/(np.sqrt(2*np.pi) * m**4 * v0**3)
    velocityfactor = (s - 4.0*(m**2)) * np.exp(-(s-4.0*m**2)/(2.0*m**2*v0**2))
    return prefactor*sigma(s)*velocityfactor

def g1 (s):
    peak_0 = 2615
    height = 3e-9
    k2 = 0.25e-3
    k1 = 0.5e-3

    ret = np.zeros(len(s),dtype=float)
    ret[s<peak_0] = height*np.exp(k1*(s[s<peak_0]-peak_0))
    ret[s>=peak_0] = height*np.exp(k2*(peak_0-s[s>=peak_0]))
    return ret


def g2 (s) :
    #height_2 / (b1 * (s - peak2)**2 + 1.0)
    return 8e-4 / (1 * (s-8317.44)**2 + 8317.44*2.5**2)

if __name__ == "__main__":
    ss = np.arange(1650,10000,1)

    pyplot.figure()
    pyplot.plot(ss,sigma(ss))
    pyplot.xlabel('Energy '+r'$s$')
    pyplot.ylabel('Cross section '+r'$\sigma$')
    pyplot.savefig('sigma.png')
    pyplot.show()

    pyplot.figure()
    pyplot.plot(ss,sigmav(ss))
    pyplot.xlabel('Energy '+r'$s$')
    pyplot.ylabel('Cross section x velocity '+r'$\sigma v$')
    #pyplot.savefig('sigmav.png')
    pyplot.plot(ss,g1(ss))
    pyplot.plot(ss,g2(ss))

    pyplot.show() 