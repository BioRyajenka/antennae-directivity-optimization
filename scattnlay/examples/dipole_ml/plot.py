#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sphereml
#dipole amplitudes
px = 1.
py = 0.
pz = 0.

wl = 455
NL = 3
RL = np.zeros(NL)
eL = np.zeros(NL+1, dtype=complex)
# Fitness = 125.09
# [0.4591807  0.37469538 0.91       0.86627232 3.07589484 3.95118756
#  2.59045949]
# RL = [0.37469538, 0.86627232,  0.91]
# eL=[3.07589484, 3.95118756, 2.59045949, 1.]
# RL = [0.37, 0.86,  0.91]
# eL=[3.08, 3.95, 2.59, 1.]
# dRd = 0.4591807*0.01 
# th=0


# Fitness = 60.765
# [2.20656848e-01 2.61279169e-02 1.53649125e-02 2.27500000e-01
#  3.44802965e+01 6.00613951e+00 4.00000000e+01]

# RL = [ 1.53649125e-02, 2.61279169e-02, 2.27500000e-01]
# eL=[3.44802965e+01, 6.00613951e+00, 4.00000000e+01, 1.]
# dRd = 2.20656848e-01*0.01 
# th=0

RL = [ 15.4, 26.1, 227.5]
eL=[3.44, 6.0, 40, 1.]
dRd = 220.656848*0.01 
th=0

wl = .455
#15.09963758950785 [0.026496548019505066, 0.020879414162895495, 0.054358920483803946, 0.05170240470232547, 19.862721369937805, 30.069986831389922, 1.0]
RL = [ 0.020879414162895495, 0.05170240470232547, 0.054358920483803946]
eL=[19.862721369937805, 30.069986831389922, 1.0, 1.]
dRd = 0.026656848*0.01
th=0


for i in range(NL): eL[i] = 2.*(2.+i) + 0.*1j
RL[0] = 0.09
for i in range(1, NL): RL[i] = RL[i-1] + 0.02
dRd = 0.001*np.sqrt(2.0);
th=np.pi

eL[NL] = 1.
data=[]
data2=[]
for i in range(138):
    Rd = dRd*i  # dipole position
    data.append([Rd,  sphereml.evaluate_directivity(RL, eL, Rd, wl, px, py, pz, th=th, ph=0., N=8)])
    data2.append([Rd,  sphereml.evaluate_directivity(RL, eL, Rd, wl, px, py, pz, th=th, ph=0., N=75)])
    # print("%8.6f %8.5f"%(Rd,  sphereml.evaluate_directivity(RL, eL, Rd, wl, px, py, pz, np.pi, 0.)))
data = np.array(data)
plt.plot(data[:,0]/wl, data[:,1], lw=2)
data2 = np.array(data2)
print(np.nanmax(data2[:,1]))
# data = np.loadtxt('results/directivity.dat')
# plt.plot(data[:,0]/wl, data[:,1], lw=0.5, color='black', marker='x')
plt.plot(data2[:,0]/wl, data2[:,1], lw=1)
for i in range(NL):
    plt.axvline(RL[i]/wl, color='red', lw=1)
plt.show()

