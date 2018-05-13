from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc
from Cython.Build import cythonize
import sys, os, shutil
import numpy

cc="icpc"
fc="ifort"
python_h_location=get_python_inc()
idir="-I./required_procedures/ -I/curc/tools/x86_64/rh6/software/gsl/1.16/intel/15.0.2/include/ -L/curc/tools/x86_64/rh6/software/gsl/1.16/intel/15.0.2/lib/ -L/curc/tools/x86_64/rh6/software/intel/15.0.2/lib/intel64/ -I/curc/tools/x86_64/rh6/software/intel/15.0.2/include/ -I/curc/tools/x86_64/rh6/software/intel/python/2.7.10/include/python2.7/ "
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
    library_dirs=["/curc/tools/x86_64/rh6/software/gsl/1.16/intel/15.0.2/lib/",
                  "/curc/tools/x86_64/rh6/software/intel/15.0.2/lib/intel64/",
                  "/curc/tools/x86_64/rh6/software/gsl/1.16/intel/15.0.2/lib/"],
    include_dirs=["./required_procedures/",
                  numpy.get_include(),
                  "/curc/tools/x86_64/rh6/software/gsl/1.16/gcc/5.1.0/include",
                  "/curc/tools/x86_64/rh6/software/intel/15.0.2/include/",
                  "./gcc-4.8.5/libstdc++-v3/include/std/"],#pointer to
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
