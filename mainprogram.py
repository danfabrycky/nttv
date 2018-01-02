# COMPILE SHARED OBJECT LIBRARY
import os
fullpath = '/home/gjgilbert/projects/koi564/markov/'
os.system('gcc -shared -Wall -o nttv.so -lgsl -lgslcblas -fPIC nttv.c')

# LOAD MODULES
import ctypes
import numpy as np
import equationsofmotion as eom
import statevectors as sv
import nttv_cwrapper

reload(eom)
reload(sv)
reload(nttv_cwrapper)


# CONVENIENCE FXNS
pi = np.pi
sin, cos, tan = np.sin, np.cos, np.tan
sqrt = np.sqrt

# VECTOR FORMATS
ps_names = ['P','T0','ecc',   'omega', 'inc','OMEGA','mp']      # pstart
po_names = ['P','T0','ecosom','esinom','inc','OMEGA','mp']      # portho
pk_names = ['a','f', 'ecc',   'omega', 'inc','OMEGA','mp']      # pkep
# state  = [x, y, z, xdot, ydot, zdot ...]                      # state vector (no mass)
# times in days; masses in M_jup; distances in AU; be careful with rad/deg for angles

# GLOBAL CONSTANTS
JOS = 9.545e-4              # M_jup/M_sol
G   = 2.9591220363e-4       # Newton's constant, in AU^3/days^2
                            # 2.9591439153e-4 (Laughlin) | 2.95912200e-4 (Dan-Eric)
# MANUAL CONTROLS
NPL   = 4
T0    = 0.
T1    = 2100.
EPOCH = 750.
MSTAR = 1.0
RSTAR = 1.0

#################
# BEGIN PROGRAM #
#################

# load and format parameter vectors
pstart = np.genfromtxt('./start.in', max_rows=NPL, skip_header=1)[:,1:]
pstart[:,[3,4,5]] *= pi/180.                                # convert angles from degrees to radians
portho = sv.start_to_ortho(pstart, 'rad')                   # switch to 'ortho' format

# calculate system mass
mpjup = portho[:,6]                                         # planet masses (M_jup)
mpsol = mpjup*JOS                                           # planet masses (M_sol)
msys  = np.array(np.cumsum(np.hstack([[MSTAR],mpsol])))     # system mass, Jacobian coordinates (M_sol)
mu    = msys/MSTAR                                          # dynamical mass (planet/star mass ratio)

# calculate keplerian and state vectors
pkep  = sv.ortho_to_kep(portho,EPOCH,msys)    
state = -sv.kep_to_state(pkep,GMtot=G*msys[1:])         # WATCH NEGATIVE SIGNS! (forward/backward integration)
pkep_ = sv.state_to_kep(-state,mpjup,GMtot=G*msys[1:])  # check self-consistency

# build tbounds vector
tbounds = np.array([T0,T1,EPOCH])

# linearize state vector
state = state.reshape(NPL*6)

# pass state vector to integrator
dummy = nttv_cwrapper.integrator(state, mu, tbounds)





















