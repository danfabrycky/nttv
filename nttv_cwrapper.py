import numpy as np
import ctypes
import os

fullpath = '/home/gjgilbert/projects/koi564/markov/'
lib = ctypes.CDLL(fullpath + 'nttv.so')


def free(pointer):
    lib.myfree.argtype = ctypes.POINTER(6*ctypes.c_double)
    
    lib.myfree(pointer)
    
    return None


def integrator(state, mu, tbounds):
    # sizes of arrays
    Nstate   = len(state)
    Nmu      = len(mu)
    Ntbounds = 3

    # set up conversion to c-readable arrays
    state_array_type    = Nstate*ctypes.c_double
    mu_array_type       = Nmu*ctypes.c_double
    tbounds_array_type  = Ntbounds*ctypes.c_double
    
    state_array         = state_array_type(*state)
    mu_array            = mu_array_type(*mu)
    tbounds_array       = tbounds_array_type(*tbounds)

    # arguments to integrator
    lib.dummy.argtypes = (ctypes.POINTER(state_array_type),    \
                               ctypes.POINTER(mu_array_type),       \
                               ctypes.POINTER(tbounds_array_type))

    # returns from integrator
    lib.dummy.restype = ctypes.POINTER(6*ctypes.c_double)
    
    integrator_ptr = lib.dummy(state_array, mu_array, tbounds_array)
    
    print integrator_ptr[0]
    
    #R = [x for x in integrator_ptr.contents]
    
    free(integrator_ptr)
        
    return 0
