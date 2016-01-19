# setup_multi.py --- setup routine for wrapping C++ H corona code

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension(
    "py_corona_multi_sim",
    ["py_corona_multi_sim.pyx"],                 # our Cython source
    language="c++",             # generate C++ code
    extra_compile_args=["-O3"],
    extra_link_args=["-lgsl","-lgslcblas","-lm","-lgfortran","-fPIC"],
    define_macros=[('SRCFNSLOC','"./source_functions/"')],
    include_dirs=["./required_procedures/"],
    extra_objects=["./required_procedures/ipbackgroundCFR_fun.o"]))
)
