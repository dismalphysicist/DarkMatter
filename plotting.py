import numpy as np
import scipy.special as sp 
import matplotlib.pyplot as plt 
import csv 

# nrows = 0 

# ss = np.zeros(nrows)
# sigma = np.zeros(nrows)
# err_sigma = np.zeros(nrows)
# sigmav = np.zeros(nrows)
# err_sigmav = np.zeros(nrows)

# with open("results.csv") as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
    
#     for row in csv_reader:
#         if nrows!=0: # not first row - titles 
#             test = 's='+row[0]+', ans='+row[1]+', err='+row[2]
#             print test 
#             ss = np.append(ss,float(row[0]))
#             sigmav = np.append(sigmav,float(row[1]))
#             err_sigmav = np.append(err_sigmav,float(row[2]))
#         nrows += 1

# print ss[-1], sigmav[-1] 

# plt.figure()
# plt.plot(ss,sigmav*0.389e9,'.')
# plt.show() 

pi = np.pi 
m = 10.
T = 0.5
v0 = np.sqrt(T/m)
print v0 

def sigma (s):
    V = A = Vtil = Atil = 1 
    M = 91.2
    gamma = 2.5

    prefactor = 1/(8*pi)
    factor1 = 1/np.sqrt(1-4*m*m/s)
    factor2 = 1/((s-M*M)*(s-M*M) + M*M * gamma*gamma)
    part1 = (V*V+A*A)*(Vtil*Vtil+Atil*Atil)*(s/3 - m*m/3)
    part2 = (V*V-A*A)*(Vtil*Vtil+Atil*Atil)*2*m*m 
    
    return prefactor*factor1*factor2*(part1+part2)

# Compare nonrelativistic and relativistic cases 

s = np.arange(4*m**2, 2000, 1)

velocity = np.sqrt(s - 4.0*(m*m))/m
prefactor = np.sqrt(2/pi)/(2*v0**3*m**3)
weight = np.sqrt(s - 4.0*(m*m)) * np.exp(-(s-4.0*m*m)/(2.0*m*m*v0*v0))

r_velocity = np.sqrt(s*(s-4*m*m))/(s-2*m*m)
r_prefactor = 1/(8*T*m**3*sp.kn(2,m/T)**2)
r_weight = np.sqrt(s - 4.0*(m*m))*(s-2*m*m)*sp.kn(1,np.sqrt(s)/T) 

sigmav = velocity*prefactor*sigma(s)*weight 
r_sigmav = r_velocity*r_prefactor*sigma(s)*r_weight 


plt.figure()
plt.plot(s, sigmav, color='b')
plt.plot(s, r_sigmav, color='r')

# plt.figure()
# plt.plot(s, velocity, color='b')
# plt.plot(s, r_velocity, color='r')

plt.show() 