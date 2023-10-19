# -*- coding: iso-8859-1 -*-

import os
import glob
from setuptools import setup
from setuptools import Extension
import versioneer

setup(name='zCluster',
      version=versioneer.get_version(),
      author='Matt Hilton + zCluster contributors',
      author_email='matt.hilton@wits.ac.za',
      packages=['zCluster'],
      package_data={'zCluster': ['data/*', 'SED/CWW/*', 'SED/BR07/*', 'SED/EAZY_v1.0/*', 
                                 'passbands/*']},
      scripts=['bin/zCluster', 'bin/zField', 'bin/zClusterBCG', 'bin/zClusterComparisonPlot'],
)
