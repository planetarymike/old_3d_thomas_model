# py_corona_sim.pyx -- Cython wrapper for H corona simulation object
# distutils: language = c++

from libcpp cimport bool
from libcpp.string cimport string

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef public api tonumpyarray(double* data, long long size) with gil:
    if not (data and size >= 0): raise ValueError
    cdef np.npy_intp dims = size
    #NOTE: it doesn't take ownership of `data`. You must free `data` yourself
    return np.PyArray_SimpleNewFromData(1, &dims, np.NPY_DOUBLE, <void*>data)

cdef extern from "nr3.h":
    cdef cppclass VecDoub:
        VecDoub() except +
        double* v
        
cdef extern from "corona_simulator.h":
    cdef cppclass corona_simulator:
        corona_simulator(bool) except +
        int nobs
        void obs_import(string, bool)
        void interp_I(double, double, VecDoub, double)
        VecDoub current_I_calc
        VecDoub I_obs
        VecDoub DI_obs
        VecDoub IPHb_model
        
cdef class Pycoronasim:
    cdef corona_simulator *thisptr #holds the reference to the cpp class
    def __cinit__(self, bool forcesim):
        self.thisptr = new corona_simulator(forcesim)
    def __dealloc__(self):
        del self.thisptr
    def obs_import(self, string obsfname, bool simulate_iph):
        self.thisptr.obs_import(obsfname, simulate_iph)
    def interp_I(self, double nH, double T, double IPHb):
        self.thisptr.interp_I(nH, T, self.thisptr.current_I_calc, IPHb)
        return tonumpyarray(self.thisptr.current_I_calc.v,self.thisptr.nobs)
    def getI_obs(self):
        return tonumpyarray(self.thisptr.I_obs.v,self.thisptr.nobs)
    def getdI_obs(self):
        return tonumpyarray(self.thisptr.DI_obs.v,self.thisptr.nobs)
    def getIPHb_model(self):	
        return tonumpyarray(self.thisptr.IPHb_model.v,self.thisptr.nobs)

