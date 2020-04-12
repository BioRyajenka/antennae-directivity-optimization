#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sphereml
import pyfde
from math import cos, pi
import random

#dipole amplitudes
px = 1.
py = 0.
pz = 0.

wl = 0.455
NL = 3
n_dim = NL*2+1
RL = np.zeros(NL)
eL = np.zeros(NL+1, dtype=complex)
eL[NL] = 1.

max_index = 30
n_it = 6000
N=50

sign="index%03.2g-N%i-NL%i-iterations%i-%05i"%(max_index, N, NL, n_it, random.randint(0,99999))

    
    
def fitness2(p):
    Rd = p[0]
    rr = np.sort(p[1:NL+1])
    if np.any(np.isclose(Rd,rr)):
        # print("problem with dipole position")
        return 0.
    for i in range(NL):
        RL[i] = rr[i]
        eL[i] = p[i+NL+1]
    D = sphereml.evaluate_directivity(RL, eL, Rd, wl, px, py, pz, th=np.pi*0., ph=0., N=N)
    if np.isnan(D): return 0.
    return D

def get_D():
    limits = [(wl*1e-3,wl*2)]
    for i in range(NL):
        limits.append((0,wl*max_ratio))
    for i in range(NL):
        limits.append((1,max_index))
    solver = pyfde.JADE(fitness2, n_dim=n_dim, n_pop=15*5, limits=limits)
    n_total = 0
    while True :
        best, fit = solver.run(n_it=20)
        n_total += 20
        print(n_total, fit)
        if n_total%1000 == 0:
            print("==> ",n_total, max_ratio, fit, best)
            with open('out2_'+sign+'.txt', 'a') as f:
                print(str(n_total)+' '+str(max_ratio)+' '+str(fit)+' '+str([i for i in best]), file=f)
        if n_total > n_it: break 
    return fit

data = []
max_ratio = 480/455
for max_ratio in np.arange(0.1,2.000005, 0.05):
    data.append([max_ratio, get_D()])
    print(data[-1])

print("--final--")
print(data)
data=np.array(data)
plt.plot(data[:,0], data[:,1])
plt.title(sign)
plt.savefig("optimization"+sign+".png")
plt.show()

# Fitness = 125.09
# [0.4591807  0.37469538 0.91       0.86627232 3.07589484 3.95118756
#  2.59045949]
