"""

Fetches all of the DECaLS catalogs needed for the given catalog and stores under DECaLSCache.

This can then be used to run zCluster on DECaLS on an HPC facility that doesn't have internet
access on the individual compute nodes (like hippo).

"""

import os
import sys
import astropy.table as atpy
from zCluster import catalogRetriever
import urllib
import time
import numpy as np

if len(sys.argv) < 2:
    print("Run: % fetchDECaLSDR8_forCatalog.py <FITS Table Path>")
else:

    tab=atpy.Table().read(sys.argv[1])
    
    # DECaLS set-up from zCluster
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
    retriever=catalogRetriever.DECaLSRetriever
    passbandSet='DECaLS'
    #zDebias=0.02    # From testing against spec-zs (don't know yet why there is a mean offset)
    maxMagError=0.2
    retrieverOptions={'maxMagError': maxMagError, 'bricksTab': bricksTab, 'DR8Tab': DR8Tab, 
                      'altCacheDir': "DECaLSCache", 'downloadOnly': True}
    
    count=0
    for row in tab:
        t0=time.time()
        count=count+1
        galaxyCatalog=catalogRetriever.DECaLSRetriever(row['RADeg'], row['decDeg'], halfBoxSizeDeg = 3.0/3600., 
                                                        optionsDict = retrieverOptions)
        t1=time.time()
        print("... %d/%d (%.3f sec) ..." % (count, len(tab), t1-t0))
        
