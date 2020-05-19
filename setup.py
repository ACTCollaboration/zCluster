# -*- coding: iso-8859-1 -*-

import os
import glob
from setuptools import setup
from setuptools import Extension
from setuptools.command.install import install
import stat
from Cython.Distutils import build_ext
import numpy
import versioneer

class OverrideInstall(install):
    """Possibly versioneer-related, but without this class, package data files are getting installed such 
    that only root can read them. This fixes the permissions so that anyone can read them.
    
    """    
    def run(self):
        mode=stat.S_IROTH
        install.run(self) 
        for filepath in self.get_outputs():
            if filepath.find("zCluster"+os.path.sep+"SED") != -1:
                os.chmod(filepath, mode)

cmdclass=versioneer.get_cmdclass()
cmdclass['build_ext']=build_ext
cmdclass['install']=OverrideInstall

setup(name='zCluster',
      version=versioneer.get_version(),
      cmdclass=cmdclass,
      url="https://github.com/ACTCollaboration/zCluster",
      author='Matt Hilton',
      author_email='matt.hilton@mykolab.com',
      classifiers=[],
      description='A code for measuring galaxy cluster photometric redshifts.',
      long_description="""A code for measuring galaxy cluster photometric redshifts. Runs on both large scale
      public survey data (e.g., SDSS) and user-supplied photometric catalogs.""",
      packages=['zCluster'],
      package_data={'zCluster': ['data/*', 'SED/CWW/*', 'SED/BR07/*', 'SED/EAZY_v1.0/*', 
                                 'passbands/*']},
      scripts=['bin/zCluster', 'bin/zField', 'bin/zClusterBCG', 'bin/zClusterComparisonPlot'],
      ext_modules=[Extension("zClusterCython", ["zCluster/zClusterCython.pyx"], include_dirs=[numpy.get_include()])],
)
