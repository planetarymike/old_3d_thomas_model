from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc
from Cython.Build import cythonize
import sys, os, shutil
import numpy

cc="icpc"
fc="ifort"
python_h_location=get_python_inc()
idir="-I./required_procedures/ -I" + os.environ['CURC_GSL_ROOT'] +"/include/ -L" + os.environ['CURC_GSL_ROOT'] +"/lib/ -L" + os.environ['CURC_INTEL_ROOT'] +"/intel64/ -I" + os.environ['CURC_GSL_ROOT'] +"/include/ -I" + os.environ['CURC_PYTHON_ROOT'] +"/include/python2.7/ "
cflags=" -lgsl -lgslcblas -lm -lfor -fPIC -lstdc++ -g -traceback"
fflags=" -lstdc++ -fPIC -nofor-main -lgsl -lgslcblas -lm -g -traceback -f77rtl -check"


c_compile = cc + " -c simulate_coronal_scan.cpp " + idir + cflags + " -O3 -o multi_corona_simulator.o"
print(c_compile)
os.system(c_compile)

f_compile = fc + " -c ./required_procedures/ipbackgroundCFR_fun.f " + idir + fflags + " -Imulti_corona_simulator.o -O3 -o cpp_corona_sim.o"
print(f_compile)
os.system(f_compile)

setup(ext_modules = cythonize(Extension(
    "py_corona_multi_sim",
    ["py_corona_multi_sim.pyx"],  # our Cython source
    language="c++",               # generate C++ code
    extra_compile_args=["-O3","-fPIC","-g","-traceback"],
    extra_link_args=["-lifcore","-lstdc++","-lgsl","-lgslcblas","-lm"],
    library_dirs=[os.environ['CURC_INTEL_ROOT'] +"/intel64/",
                  os.environ['CURC_GSL_ROOT'] +"/lib"],
    include_dirs=["./required_procedures/",
                  numpy.get_include(),
                  os.environ['CURC_GSL_ROOT'] +"/include",
                  os.environ['CURC_GSL_ROOT'] +"/include"],
#                  "./gcc-4.8.5/libstdc++-v3/include/std/"],#pointer to
                                                           #a manually
                                                           #installed
                                                           #version of
                                                           #GCC
                                                           #headers,
                                                           #the intel
                                                           #compiler
                                                           #is
                                                           #incompatible
                                                           #with
                                                           #versions >
                                                           #4.8
    extra_objects=["./cpp_corona_sim.o"]))
)
