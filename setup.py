# -*- coding: iso-8859-1 -*-

import os
import glob
from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
#import popen2
#import numpy

setup(name='zCluster',
      version="git",
      url=None,
      author='Matt Hilton',
      author_email='matt.hilton@mykolab.com',
      classifiers=[],
      description='A code for measuring galaxy cluster photometric redshifts.',
      long_description="""A code for measuring galaxy cluster photometric redshifts. Runs on both large scale
      public survey data (e.g., SDSS) and user-supplied photometric catalogs.""",
      packages=['zCluster'],
      package_data={'zCluster': ['data/*', 'SED/CWW/*', 'SED/BR07/*', 'SED/EAZY_v1.0/*', 
                                 'SED/PEGASE2.0/*', 'SED/SWIRETemplateLibrary/*', 
                                 'passbands/*']},
      scripts=['bin/zCluster', 'bin/zClusterBCG', 'bin/zClusterComparisonPlot'],
      #cmdclass={'build_ext': build_ext},
      #ext_modules=[Extension("filtersCython", ["zCluster/filtersCython.pyx"], extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp'], include_dirs=[numpy.get_include()])]

)
