from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
import distutils.sysconfig as ds

setup(name='PyBLS',
      version='0.1',
      description='Python wrapper for the BLS routine.',
      author='Hannu Parviainen',
      author_email='hpparvi@gmail.com',
      url='',
      extra_options = ['-fopenmp'],
      package_dir={'pybls':'src'},
      packages=['pybls'],
      ext_modules=[Extension('pybls.blsf', ['src/bls.f90'], libraries=['gomp','m'])]
     )
