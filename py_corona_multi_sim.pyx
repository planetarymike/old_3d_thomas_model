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
        
cdef extern from "multi_corona_simulator.h":
    cdef cppclass corona_observation:
        corona_observation(string, bool) except +
        int nobs
        VecDoub current_I_calc
        VecDoub I_obs
        VecDoub DI_obs
        VecDoub IPHb_model        
    cdef cppclass corona_simulator:
        corona_simulator(bool,bool) except +
        int nobsdata
        corona_observation* allobsdata
        void interp_I(int, double, double, double)
        void load_obs(string)
        void add_obs(string)
        
cdef class Pycorona_multi_sim:
    cdef corona_simulator *thisptr #holds the reference to the cpp class
    def __cinit__(self, bool forcesim, bool iphsim):
        self.thisptr = new corona_simulator(forcesim,iphsim)
    def __dealloc__(self):
        del self.thisptr
    def obs_import(self, string obsfname):
        self.thisptr.add_obs(obsfname)
    def interp_I(self, int idata, double nH, double T, double IPHb):
        self.thisptr.interp_I(idata, nH, T, IPHb)
        return tonumpyarray(self.thisptr.allobsdata[idata].current_I_calc.v,self.thisptr.allobsdata[idata].nobs)
    def getI_obs(self, int idata):
        return tonumpyarray(self.thisptr.allobsdata[idata].I_obs.v,self.thisptr.allobsdata[idata].nobs)
    def getdI_obs(self,int idata):
        return tonumpyarray(self.thisptr.allobsdata[idata].DI_obs.v,self.thisptr.allobsdata[idata].nobs)
    def getIPHb_model(self,int idata):	
        return tonumpyarray(self.thisptr.allobsdata[idata].IPHb_model.v,self.thisptr.allobsdata[idata].nobs)
    def getnobsdata(self):
        return self.thisptr.nobsdata

