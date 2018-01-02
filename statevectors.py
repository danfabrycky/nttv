###
#   Functions for converting between different formats for state vectors
#
#   'start':  p = [P, T0, ecc, inc, OMEGA, omega, mp ...]
#   'ortho':  p = [P, T0, ecosom, esinom, inc, OMEGA, mp ...]
#   'kepler': p = [a, ecc, inc, omega, OMEGA, f, mp ...]
#   'state':  p = [x, y, z, xdot, ydot, zdot ...]
###

import numpy as np
import equationsofmotion as eom

reload(eom)

pi = np.pi
sin, cos, tan = np.sin, np.cos, np.tan
sqrt = np.sqrt

###

def start_to_ortho(pstart, angles):
    '''
    convert from 'start' to 'ortho' format or vice versa
    
    pstart  = parameters, (NPL x 7) array
    angles  = can be 'deg' or 'rad'
    '''
    # set angle conversion factor for degrees/radians
    if angles == 'deg': ac = pi/180.
    if angles == 'rad': ac = 1.0

    # pull values from parameter array
    ps = np.swapaxes(pstart,0,1)
    
    P, T0, ecc, omega, inc, OMEGA, mp = np.copy(ps)
    esinom = ecc*sin(omega*ac)
    ecosom = ecc*cos(omega*ac)

    return np.swapaxes(np.array([P, T0, ecosom, esinom, inc, OMEGA, mp]),0,1)



def ortho_to_start(portho, angles):
    '''
    convert from 'start' to 'ortho' format or vice versa
    
    portho  = parameters, (NPL x 7) array
    angles  = can be 'deg' or 'rad'
    '''
    # set angle conversion factor for degrees/radians
    if angles == 'deg': ac = pi/180.
    if angles == 'rad': ac = 1.0

    # pull values from parameter array
    ps = np.swapaxes(portho,0,1)
    
    P, T0, ecosom, esinom, inc, OMEGA, mp = np.copy(ps)
    ecc   = sqrt(ecosom**2 + esinom**2)
    omega = np.arctan2(esinom,ecosom)/ac

    return np.swapaxes(np.array([P, T0, ecc, omega, inc, OMEGA, mp]),0,1)



def ortho_to_kep(portho, epoch, msys):
    '''
    convert from 'ortho' format to 'kepler' format

    portho  = parameters, (NPL x 7) array
    epoch   = transit epoch (days)
    msys    = mass of system (M_sol; Jacobian coordinates), length(NPL)
    
    input angles must be in RADIANS

    '''
    # set constants
    ghere = 2.9591220363e-4             # Newtons constant

    # pull values from parameter array
    ps = np.swapaxes(portho,0,1)
    
    P, T0, ecosom, esinom, inc, OMEGA, mp = np.copy(ps)
    ecc   = sqrt(ecosom**2 + esinom**2)
    omega = np.arctan2(esinom,ecosom)
        
    pomega = OMEGA + omega
    lamb0  = eom.getlamb(-omega+pi/2., ecc, pomega)                     # mean longitude
    M0     = lamb0 - pomega                                             # mean anomaly
    Me     = M0 + 2.*pi*(epoch-T0)/P                                    # mean anomaly at epoch

    lambe  = (Me+pomega) % (2*pi)                                       # mean longitude at epoch
    lambe[lambe>pi]  -= 2*pi                                            # range -pi < lambe < +pi
    lambe[lambe<-pi] += 2*pi
 
    f      = eom.getf(lambe, ecc, pomega)                               # true anomaly
    a      = np.power(ghere*msys[1:],1./3) * np.power(P/(2.*pi),2./3)   # semimajor axis
    
    return np.swapaxes(np.array([a, f, ecc, omega, inc, OMEGA, mp]),0,1)



def rotatekep(x,y,inc,omega,OMEGA):
    cosom = cos(omega); sinom = sin(omega)
    cosOM = cos(OMEGA); sinOM = sin(OMEGA)
    cosi  = cos(inc);   sini  = sin(inc)

    x1 = cosom*x - sinom*y
    y1 = sinom*x + cosom*y
    z1 = np.zeros_like(x1)

    x2 = x1
    y2 = cosi*y1 - sini*z1
    z2 = sini*y1 + cosi*z1

    x3 = cosOM*x2 - sinOM*y2
    y3 = sinOM*x2 + cosOM*y2
    z3 = z2
        
    return np.array([x3,y3,z3])



def kep_to_state(pkep,GMtot=1.0):
    '''
    Inputs: pkep  = Keplerian parameter vector [a, f, ecc, inc, omega, OMEGA, mp]
            GMtot = system total mass
            
    Outputs: state = position & velocity state vector [x,y,z,xdot,ydot,zdot]
    '''
    pkeps = np.swapaxes(pkep,0,1)

    a, f, ecc, omega, inc, OMEGA, mp = pkeps
    
    cosf  = cos(f)
    sinf  = sin(f)    
    a1me2 = a*(1.0-ecc**2)
    
    r = a1me2/(1.0+ecc*cosf)
       
    x0 = r*cosf
    y0 = r*sinf
    
    vx0 = -sqrt(GMtot/a1me2)*sinf
    vy0 =  sqrt(GMtot/a1me2)*(ecc+cosf)

    xyz = np.swapaxes(rotatekep(x0,y0,inc,omega,OMEGA),0,1)
    xyzdot = np.swapaxes(rotatekep(vx0,vy0,inc,omega,OMEGA),0,1)
    
    return np.hstack([xyz, xyzdot])


def state_to_kep(state, mp, GMtot=1.0):
    '''
    Inputs: state = position & velocity state vector [x,y,z,xdot,ydot,zdot]
            GMtot = system total mass
            
    Outputs: pkep  = Keplerian parameter vector [a, f, ecc, inc, omega, OMEGA, mp]
    
    For all equations see Murray & Dermott chapter 2.8
    '''
    x,y,z,vx,vy,vz = np.swapaxes(state,0,1)
    
    rsq = x**2 + y**2 + z**2
    r   = sqrt(rsq)
    vsq = vx**2 + vy**2 + vz**2
    v   = sqrt(vsq)    
    hx  = y*vz - z*vy
    hy  = z*vx - x*vz
    hz  = x*vy - y*vx
    hsq = hx**2 + hy**2 + hz**2
    h   = sqrt(hsq)

    rdot = sqrt(abs(vsq - hsq/rsq))
    rdot[vsq <= (hsq/rsq)] = 0.0

    rrdot = x*vx + y*vy + z*vz
    rdot[rrdot <= 0] *= -1
    
    a   = (2./sqrt(rsq) - vsq/GMtot)**-1
    ecc = sqrt(1. - hsq/(GMtot*a))
    inc = np.arccos(hz/h)
    
    sini  = sin(inc)
    
    sinOM = hx/(h*sini)
    cosOM = -hy/(h*sini)
    OMEGA = np.arctan2(sinOM,cosOM)
    
    sinomf =  z/(r*sini)
    cosomf = (1./cos(OMEGA))*(x/r + sin(OMEGA)*sinomf*hz/h)
    omf    = np.arctan2(sinomf,cosomf)

    OMEGA[sini==0] = 0.0
    omf[sini==0]   = np.arctan2(y/r,x/r)[sini==0]
    
    sinf = a*(1.-ecc**2)/(h*ecc)*rdot
    cosf = 1./ecc *(a*(1.-ecc**2)/r - 1)
    f    = np.arctan2(sinf,cosf)
    
    omega = (omf - f) % (2*pi)      # -pi < omega < pi
    omega[omega>pi]  -= (2*pi)
    omega[omega<-pi] += (2*pi)
    
    omega[ecc==0] = 0.0
    f[ecc==0] = omf[ecc==0]
    
    return np.swapaxes(np.array([a, f, ecc, omega, inc, OMEGA, mp]),0,1)
