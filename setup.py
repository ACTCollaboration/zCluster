# -*- coding: iso-8859-1 -*-

import os
import glob
from setuptools import setup
from setuptools import Extension
from setuptools.command.install import install
import stat
#from Cython.Distutils import build_ext
import numpy
import versioneer

cmdclass=versioneer.get_cmdclass()
#cmdclass['build_ext']=build_ext

setup(name='zCluster',
      version=versioneer.get_version(),
      cmdclass=cmdclass,
      url="https://github.com/ACTCollaboration/zCluster",
      author='Matt Hilton + zCluster contributors',
      author_email='matt.hilton@mykolab.com',
      classifiers=[],
      description='A code for measuring galaxy cluster photometric redshifts.',
      long_description="""A code for measuring galaxy cluster photometric redshifts. Runs on both large scale
      public survey data (e.g., SDSS) and user-supplied photometric catalogs.""",
      packages=['zCluster'],
      package_data={'zCluster': ['data/*', 'SED/CWW/*', 'SED/BR07/*', 'SED/EAZY_v1.0/*', 
                                 'passbands/*']},
      scripts=['bin/zCluster', 'bin/zField', 'bin/zClusterBCG', 'bin/zClusterComparisonPlot'],
      #ext_modules=[Extension("zClusterCython", ["zCluster/zClusterCython.pyx"], include_dirs=[numpy.get_include()])],
      install_requires=["astropy >= 4.0",
                        "numpy >= 1.19",
                        "matplotlib >= 2.0",
                        "astLib >= 0.11.7",
                        "scipy >= 1.0",
                        #"cython",
                        "requests"]
)
