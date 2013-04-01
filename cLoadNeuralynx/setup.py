'''
To compile Cython run:
 python setup.py build_ext --inplace
'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np        # Added to avoid "arrayobject.h" error

sourcefiles = ['cLoadNeuralynx.pyx']


ext_modules = [Extension(
    name="cLoadNeuralynx",
    sources=["cLoadNeuralynx.pyx"],
    extra_objects=["loadNeuralynx.o"],
    )]

setup(
    name = 'cLoadNeuralynx',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],  # Added to avoid "arrayobject.h" error
    ext_modules = ext_modules
    )
