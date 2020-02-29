import sphereml
import numpy as np

def evaluate(ind, layers=2):
    assert len(ind)==layers*3+1
    
    r = ind[-1]
    x = ind[:layers]
    m_ = ind[layers:3*layers][::2]
    m_j = ind[layers:3*layers][1::2]
    m = map(lambda t: t[0]+t[1]*1j, zip(m_, m_j))
    
    m = [z for _,z in sorted(zip(x,m), key=lambda pair: pair[0])]
    x = sorted(x)
    
    m = np.array(m)
    x = np.array(x)
    
    px = 1.
    py = 0.
    pz = 0.

    wl = 0.455
    NL = layers
    RL = np.zeros(NL)
    eL = np.zeros(NL+1, dtype=complex)
    eL[:NL] = m[:]
    eL[NL] = 1.
    RL[:] = x[:] / 1000.0
    Rd = r / 1000.0
    
    res = sphereml.evaluate_directivity(RL, eL, Rd, wl, px, py, pz, np.pi, 0.)
    if np.isnan(res): 
        print("!!! ERROR {}".format(ind))
        print(" eL:{}".format(eL))
        print(" RL:{}".format(RL))
        print("  x:{}".format(x))
    
    return res