import numpy as np
import ctypes
import os


# compile and load library
fullpath = '/home/gjgilbert/cwrapper/tutorial_ode/'
os.system('gcc -shared -Wall -o vdp.so -lgsl -lgslcblas -fPIC vdp.c')

lib = ctypes.CDLL('/home/gjgilbert/cwrapper/tutorial_ode/vdp.so')


def free(pointer):
    lib.myfree.argtype = ctypes.POINTER(3*ctypes.c_double)
    
    lib.myfree(pointer)
    
    return None


def solver():
    lib.solver.restype = ctypes.POINTER(6*ctypes.c_double)
    
    solver_ptr = lib.solver()
    
    R = [x for x in solver_ptr.contents]
    
    free(solver_ptr)
    
    return R

s = solver()

print s
    
    
