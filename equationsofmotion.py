# Equations from Murray & Dermott (1999)
# angle inputs to all functions **MUST** be in radians

import numpy as np

pi = np.pi
sin, cos, tan = np.sin, np.cos, np.tan
arcsin, arccos, arctan = np.arcsin, np.arccos, np.arctan
sqrt = np.sqrt


def getE(f,ecc):
    return 2.*arctan(sqrt((1.-ecc)/(1.+ecc))*tan(f/2.))     # 2.46
    

def getlamb(f,ecc,pomega):
    E = getE(f,ecc)                     # 2.46
    lamb = pomega + E - ecc*sin(E)      # 2.52 & 2.53
    return lamb
    

def getf(lamb,ecc,pomega):
    M = lamb - pomega                   # 2.53
    E = M                               # 2.64 for slight improvement when using auto-convergence
    
    for i in range(25):
        E = M + ecc*sin(E)              # 2.52, Newton-Raphson would be faster
        
    return 2.*arctan( sqrt(1.+ecc)/sqrt(1.-ecc) * tan(E/2.) )

