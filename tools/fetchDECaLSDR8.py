"""

Fetches all of the DECaLS catalogs and stores under DECaLSCache.

This can then be used to run zCluster on DECaLS on an HPC facility that doesn't have internet
access on the individual compute nodes (like hippo).

This should clock in at ~2.6 TB, so probably you want to use fetchDECaLSDR8_forCatalog.py instead...

"""

import os
import sys
import astropy.table as atpy
import urllib
import numpy as np

bricksCacheDir="DECaLSCache"
os.makedirs(bricksCacheDir, exist_ok = True)
bricksPath=bricksCacheDir+os.path.sep+"survey-bricks.fits.gz"
if os.path.exists(bricksPath) == False:
    print("... fetching and caching DECaLS survey-bricks.fits.gz ...")
    urllib.request.urlretrieve("http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz", bricksPath)
bricksDR8Path=bricksCacheDir+os.path.sep+"survey-bricks-dr8-south.fits.gz"
if os.path.exists(bricksDR8Path) == False:
    print("... fetching and caching DECaLS survey-bricks-dr8-south.fits.gz ...")
    urllib.request.urlretrieve("https://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/south/survey-bricks-dr8-south.fits.gz", bricksDR8Path)
bricksTab=atpy.Table().read(bricksPath)
DR8Tab=atpy.Table().read(bricksDR8Path)
DR8Tab.rename_column("brickname", "BRICKNAME")

basePath="https://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/south/tractor/"
matchTab=atpy.join(bricksTab, DR8Tab, keys = 'BRICKNAME')
cacheDir=bricksCacheDir
count=0
for mrow in matchTab:
    count=count+1
    print("... %d/%d ..." % (count, len(matchTab)))
    url=basePath+"%03d" % np.floor(mrow['RA'])
    url=url+os.path.sep+"tractor-"+mrow['BRICKNAME']+".fits"
    fileName=cacheDir+os.path.sep+"tractor-%s.fits" % (mrow['BRICKNAME'])
    if os.path.exists(fileName) == False:
        print("... retrieving tractor catalog from web ...")
        urllib.request.urlretrieve(url, filename = fileName)
