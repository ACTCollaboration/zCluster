"""

Fetches all of the DECaLS catalogs needed for the given catalog and stores under a user-specified
cache directory.

"""

import os
import sys
import astropy.table as atpy
from zCluster import retrievers
import urllib
import time
import numpy as np

if len(sys.argv) < 5:
    print("Run: % fetchDECaLSForCatalog.py <catalog> <DR8|DR9> <path to cache> <half box size arcmin>")
else:

    tab=atpy.Table().read(sys.argv[1])
    DR=sys.argv[2]
    if DR not in ['DR8', 'DR9']:
        raise Exception("Choose either 'DR8' or 'DR9'")
    cacheDir=sys.argv[3]
    halfBoxSizeDeg=float(sys.argv[4])/60
    
    os.makedirs(cacheDir, exist_ok = True)
    database="DECaLS%s" % (DR)
    retriever, retrieverOptions, passbandSet=retrievers.getRetriever(database, maxMagError = 0.2)
    
    count=0
    for row in tab:
        t0=time.time()
        count=count+1
        galaxyCatalog=retriever(row['RADeg'], row['decDeg'], halfBoxSizeDeg = halfBoxSizeDeg, 
                                optionsDict = retrieverOptions)
        t1=time.time()
        print("... %d/%d (%.3f sec) ..." % (count, len(tab), t1-t0))
        
