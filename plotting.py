import numpy as np
import matplotlib.pyplot as plt 
import csv 

# Open file 
nrows = 0 

ss = np.zeros(nrows)
sigma = np.zeros(nrows)
err_sigma = np.zeros(nrows)
sigmav = np.zeros(nrows)
err_sigmav = np.zeros(nrows)

with open("results.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    
    for row in csv_reader:
        if nrows!=0: # not first row - titles 
            test = 's='+row[0]+', ans='+row[1]+', err='+row[2]
            print test 
            ss = np.append(ss,float(row[0]))
            sigmav = np.append(sigmav,float(row[1]))
            err_sigmav = np.append(err_sigmav,float(row[2]))
        nrows += 1

print ss[-1], sigmav[-1] 

plt.figure()
plt.plot(ss,sigmav*0.389e9,'.')
plt.show() 