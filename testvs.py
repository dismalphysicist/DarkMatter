import numpy as np
import matplotlib.pyplot as pyplot 


def abs2(arr):
    if arr.size==3:
        return sum(arr*arr)
    if arr.size==4:
        return arr[0]**2 - arr[1]**2 - arr[2]**2 - arr[3]**2 

m = 10

vel1 = np.array([0.,0.,1.])/3e2    # Defined without gamma 
gamma1 = 1/np.sqrt(1-abs2(vel1))
mom1 = np.array([m*gamma1,m*gamma1*vel1[0],m*gamma1*vel1[1],m*gamma1*vel1[2]])

vel2 = np.array([0.,0.,-0.5])/3e2   # Defined without gamma 
gamma2 = 1/np.sqrt(1-abs2(vel2))
mom2 = np.array([m*gamma2,m*gamma2*vel2[0],m*gamma2*vel2[1],m*gamma2*vel2[2]])

s = abs2(mom1+mom2)     # Definition of s 

v_CM = np.sqrt(abs2(vel1+vel2))

v1_CM = (np.sqrt(abs2(vel1))-v_CM)/(1-np.sqrt(abs2(vel1))*v_CM)

vrel = np.sqrt(abs2(vel1-vel2))     # Lab frame 
#vrel_CM = vrel*(1-v_CM**2)/(1+v_CM**2*(np.sqrt(abs2(vel1))*np.sqrt(abs2(vel2))-0.5))   # CoM frame, if v_CM and v_rel are parallel 
vrel_CM = (v1_CM**2 - v_CM**2)/(1+v1_CM*v_CM)**2

sp = 4*m**2*(1 + vrel**2/(4-vrel**2))   # Is this valid generally? 
spp = 4*m**2 + m**2*vrel**2   # Non-relativistic limit of sp 

# sp_ = 4*m**2*(1 + v1CM2/(1-v1CM2))
# spp_ = 4*m**2 + 2*m**2*v1CM2    # Doesn't work

s_ = 4*m**2*(1 + vrel_CM**2/(4-vrel_CM**2))     # This should be valid generally 


print "\t\t\ts = ", s
print "\nIn terms of v_rel (lab),s = ", sp 
print "Non-relativistic limit: s = ", spp 
print "\nIn terms of v_rel (CM), \ts = ", s_ 