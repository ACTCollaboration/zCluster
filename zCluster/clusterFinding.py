"""
    Routines for finding cluster redshift from p(z) distributions of galaxies in catalog

    ---
    
    Copyright 2017 Matt Hilton (matt.hilton@mykolab.com)
    
    This file is part of zCluster.

    zCluster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    zCluster is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with zCluster.  If not, see <http://www.gnu.org/licenses/>.
    
"""

import os
import sys
import numpy as np
from astLib import *
from scipy import interpolate
from scipy import stats
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from PIL import Image
from PIL import ImageDraw
import nemoCython
import pylab as plt
import IPython

#------------------------------------------------------------------------------------------------------------
def getPixelAreaDeg2Map(mapData, wcs):
    """Returns a map of pixel area in square degrees.
    
    """
    
    # Get pixel size as function of position
    pixAreasDeg2=[]
    RACentre, decCentre=wcs.getCentreWCSCoords()
    x0, y0=wcs.wcs2pix(RACentre, decCentre)
    x1=x0+1
    for y0 in range(mapData.shape[0]):
        y1=y0+1
        ra0, dec0=wcs.pix2wcs(x0, y0)
        ra1, dec1=wcs.pix2wcs(x1, y1)
        xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
        yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
        pixAreasDeg2.append(xPixScale*yPixScale)
    pixAreasDeg2=np.array(pixAreasDeg2)
    pixAreasDeg2Map=np.array([pixAreasDeg2]*mapData.shape[1]).transpose()
    
    return pixAreasDeg2Map   

#-------------------------------------------------------------------------------------------------------------
def makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, weightsType, minDistanceMpc = 0.0, maxDistanceMpc = 1.0, 
                   applySanityCheckRadius = True, sanityCheckRadiusArcmin = 1.0, areaMask = None, wcs = None):
    """Make a N(z) distribution in the direction of the postion (usually of a cluster) given by RADeg, decDeg.
    This is constructed from the p(z) distributions of all galaxies in the catalog, with radial weights applied
    as per weightsType. So, for 'flat1Mpc', the radial weight is 1 if within a projected distance of 1 Mpc of 
    RADeg, decDeg, and 0 otherwise.
    
    zPriorMax can be used to reign in any spuriously high redshifts given our prior knowledge of the depth of a
    given photometric survey. However, this is now largely redundant as the mag-based prior does most of the
    work (see PhotoRedshiftEngine.py).
    
    NOTE: minDistanceMpc only applies to 'flat' weightsType (for local background estimate)
    
    """

    # For speed, we make a big array of all p(z) here
    pzArray=[]
    for obj in catalog:
        if 'pz' in list(obj.keys()):
            pzArray.append(obj['pz'])
    pzArray=np.array(pzArray)
    
    # Make big array of angular separation for each galaxy
    # Get EAZY odds too, in case useful for weighting
    gOdds, gRedshifts, angSepArray, tanArray=extractArraysFromGalaxyCatalog(catalog, RADeg, decDeg)
    
    # Sanity check: if we didn't find a decent number of galaxies close to the cluster centre, 
    # the cluster is probably not actually in the optical catalog footprint (e.g., just outside edge of S82)
    if applySanityCheckRadius == True and np.sum(np.less(np.degrees(tanArray), sanityCheckRadiusArcmin/60.0)) == 0:
        print("... no galaxies within %.1f' of cluster position ..." % (sanityCheckRadiusArcmin))
        return None

    # Assume same z grid for all objects, calculate DAs
    for obj in catalog:
        if 'pz_z' in list(obj.keys()):
            zArray=obj['pz_z']
            break
    DAArray=[]
    for z in zArray:
        DAArray.append(astCalc.da(z))
    DAArray=np.array(DAArray)
        
    # NOTE: radial weights should have max value 1, min value 0
    # Koester, Postman et al. style radial weights (projected 2d NFW profile)
    if weightsType == 'NFW':
        rs=0.15 # rs=R200/c c=concentration - Koester uses 150 kpc
        r=np.linspace(0.0, 10.0, 2000)
        x=r/rs
        fx=np.zeros(x.shape)
        mask=np.greater(x, 1)
        fx[mask]=1-(2.0/np.sqrt(x[mask]**2-1))*np.arctan(np.sqrt((x[mask]-1)/(x[mask]+1)))
        mask=np.less(x, 1)
        fx[mask]=1-(2.0/np.sqrt(1-x[mask]**2))*np.arctanh(np.sqrt((1-x[mask])/(x[mask]+1)))
        mask=np.equal(x, 1)
        fx[mask]=0
        mask=np.greater(x, 20)
        fx[mask]=0
        sigmax=np.zeros(r.shape)
        mask=np.greater(r, rs)
        sigmax[mask]=(2*rs*fx[mask])/(x[mask]**2-1)
        mask=np.less(r, rs)
        sigmax[mask]=sigmax.max()
        sigmax=sigmax/sigmax.max()
        tckSigmax=interpolate.splrep(r, sigmax)
        rMpc=(np.array([tanArray]*DAArray.shape[0]).transpose())*DAArray
        rWeights=[]
        for row in rMpc:
            rWeights.append(interpolate.splev(row, tckSigmax))
        rWeights=np.array(rWeights)
        rWeights[np.less(rWeights, 0)]=0.0
        # We have to match area to background: so, zero weight all galaxies outside 1 Mpc radius
        rWeights[np.greater(rMpc, maxDistanceMpc)]=0.0
    
    elif weightsType == 'flat':
    
        # Calculate weights - set to zero for objects outside of fiducial projected distance
        # 1 inside 1 Mpc radius, 0 otherwise
        # Use this for background (the 1 Mpc radius just normalises area
        rMpc=(np.array([tanArray]*DAArray.shape[0]).transpose())*DAArray
        rWeights=np.ones(rMpc.shape)
        rWeights[np.greater(rMpc, maxDistanceMpc)]=0.0
        rWeights[np.less(rMpc, minDistanceMpc)]=0.0
        
    elif weightsType == 'radial':
        
        # Calculate weights - set to zero for objects outside of fiducial projected distance
        # 1 if at centre coords, falls linearly with increasing radius
        rMpc=(np.array([tanArray]*DAArray.shape[0]).transpose())*DAArray
        rWeights=1.0-rMpc/maxDistanceMpc
        rWeights[np.greater(rMpc, maxDistanceMpc)]=0.0
    
    else:
        raise Exception("didn't understand weightsType")
    
    # Simplest odds weighting scheme is to cut everything with odds < some value
    # Whatever happens, we definitely want to cut anything with odds == 0, as its z will default to zMin
    # NOTE: leave set at gOddsCut=0.5 - this fixed null test problems and didn't negatively impact z accuracy
    gOddsCut=0.5
    gOddsWeights=np.ones(gOdds.shape)
    mask=np.less(gOdds, gOddsCut)
    gOddsWeights[mask]=0.0    
    gOddsWeights=[gOddsWeights]*rWeights.shape[1]
    gOddsWeights=np.array(gOddsWeights).transpose()
    #gOddsWeights[:]=1 # uncomment to disable gOdds weights
    
    # NOTE: zPrior applied in calculating cluster redshift later, not in counting galaxies here
    Nz_r=rWeights.sum(axis = 0)                                             # Number of galaxies within 1 Mpc at each z (IF all weights were 1)
    Nz_r_odds_total=np.greater(rWeights*gOddsWeights, 0).sum(axis = 0)      # Actual number of galaxies within 1Mpc at each z regardless of weight (if using NFW, radial - for flat, Nz_total == Nz_r)
    Nz_r_odds=np.sum(rWeights*gOddsWeights, axis = 0)   # Number of galaxies within 1 Mpc at each z AND odds > cut
    NzWeightedSum=np.sum(pzArray*rWeights*gOddsWeights, axis = 0)
    
    # We may as well keep track of area as well
    if areaMask is None:
        areaMpc2=np.array([np.pi*maxDistanceMpc**2 - np.pi*minDistanceMpc**2]*len(zArray))
    else:
        areaMap=getPixelAreaDeg2Map(areaMask, wcs)
        rDegMap=np.zeros(areaMask.shape)
        rDegMap, xRange, yRange=nemoCython.makeDegreesDistanceMap(rDegMap, wcs, RADeg, decDeg, 2.0)
        rDegMap=rDegMap[yRange[0]:yRange[1], xRange[0]:xRange[1]]
        areaMap=areaMask[yRange[0]:yRange[1], xRange[0]:xRange[1]]*areaMap[yRange[0]:yRange[1], xRange[0]:xRange[1]]
        areaMpc2=[]
        count=0
        for DA in DAArray:
            areaScaling=1/np.power(np.degrees(np.arctan(1.0/DA)), 2)
            max_rDegCut=np.degrees(np.arctan(maxDistanceMpc/DA))
            min_rDegCut=np.degrees(np.arctan(minDistanceMpc/DA))
            areaDeg2=areaMap[np.logical_and(rDegMap >= min_rDegCut, rDegMap < max_rDegCut)].sum()
            areaMpc2.append(areaScaling*areaDeg2)
            count=count+1
        areaMpc2=np.array(areaMpc2)

    # If we want p(z), we can divide pzWeightedSum by Nz_r_odds
    return {'NzWeightedSum': NzWeightedSum, 'Nz': Nz_r_odds, 'Nz_total': Nz_r_odds_total, 'zArray': zArray, 
            'DAArray': DAArray, 'areaMpc2': areaMpc2}
    
#-------------------------------------------------------------------------------------------------------------
def estimateClusterRedshift(RADeg, decDeg, catalog, zPriorMin, zPriorMax, weightsType, maxRMpc, 
                            zMethod, maskMap = None, maskWCS = None, sanityCheckRadiusArcmin = 1.0, 
                            bckCatalog = [], bckAreaDeg2 = None, zDebias = None):
    """This does the actual work of estimating cluster photo-z from catalog.
    
    Assumes each object has keys 'pz' (p(z), probability distribution), 'pz_z' (corresponding redshifts at 
    each point in p(z)
    
    If bckCatalog and bckAreaDeg2 are given, then these are used for the background estimation. If not, then
    a projected 3-4 Mpc annulus is used.
    
    If zDebias is not None, the final output redshift ('z') will be:
        
        z_final = z_initial + zDebias*(1+z_initial) 
    
    This is obviously a fudge, so use with caution (only used for some surveys - see bin/zCluster).
            
    """
    
    print("... estimating cluster photo-z ...")
    
    # Initial sanity check: if we didn't find a decent number of galaxies close to the cluster centre, 
    # the cluster is probably not actually in the optical catalog footprint (e.g., just outside edge of S82)
    gOdds, gRedshifts, angSepArray, tanArray=extractArraysFromGalaxyCatalog(catalog, RADeg, decDeg)
    if np.sum(np.less(np.degrees(tanArray), sanityCheckRadiusArcmin/60.0)) == 0:
        print("... no galaxies within %.1f' of cluster position ..." % (sanityCheckRadiusArcmin))
        return None

    # Automated area mask calculation, using only catalogs
    # Probably not super accurate but better than nothing - will at least take care of survey boundaries
    if maskMap is None and maskWCS is None:
        areaMask, wcs=estimateAreaMask(RADeg, decDeg, catalog)
    else:
        areaMask, wcs=maskMap, maskWCS
    
    # Here, we do a weighted average of the p(z) distribution
    # The weighting is done according to distance from the cluster position at the z corresponding to a
    # given p(z). 
    clusterNzDict=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, weightsType, maxDistanceMpc = 1.0,
                                 areaMask = areaMask, wcs = wcs)
    zArray=clusterNzDict['zArray']
    if clusterNzDict['NzWeightedSum'].sum() == 0:
        print("... no galaxies in n(z) - skipping...")
        return None
    
    # This is how we have been calculating redshifts
    # We only need normFactor to make sure odds is correct
    pzWeightedMean=clusterNzDict['NzWeightedSum']/clusterNzDict['Nz']
    pzWeightedMean[np.isnan(pzWeightedMean)]=0.0
    normFactor=np.trapz(pzWeightedMean, clusterNzDict['zArray'])
    pzWeightedMean=pzWeightedMean/normFactor
    try:
        z, odds, zOdds=calculateRedshiftAndOdds(pzWeightedMean, clusterNzDict['zArray'], dzOdds = 0.05, method = zMethod, zPriorMax = zPriorMax, zPriorMin = zPriorMin)
    except:
        print("Hmm - failed in z, odds calc: what happened?")
        IPython.embed()
        sys.exit()
        
    # Background - delta (density contrast) measurement
    clusterNzDictForSNR=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, 'flat', maxDistanceMpc = maxRMpc)
    zIndex=np.where(zArray == z)[0][0]
    if len(bckCatalog) > 0:
        # Global
        bckNzDict=makeWeightedNz(RADeg, decDeg, bckCatalog, zPriorMax, 'flat', minDistanceMpc = 3.0, maxDistanceMpc = 1e9, 
                                 applySanityCheckRadius = False)
        bckNzDict['areaMpc2']=((np.radians(np.sqrt(bckAreaDeg2))*bckNzDict['DAArray'])**2)
    else:
        # Local background - 0.6 deg radius gives us out to 4 Mpc radius at z = 0.1
        bckNzDict=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, 'flat', minDistanceMpc = 3.0, maxDistanceMpc = 4.0,
                                 areaMask = areaMask, wcs = wcs)
    # Both
    bckAreaNorm=clusterNzDictForSNR['areaMpc2'][zIndex]/bckNzDict['areaMpc2'][zIndex]
    bckSubtractedCount=clusterNzDictForSNR['NzWeightedSum'][zIndex]-bckAreaNorm*bckNzDict['NzWeightedSum'][zIndex]
    delta=bckSubtractedCount/(bckAreaNorm*bckNzDict['NzWeightedSum'][zIndex])
    errDelta=np.sqrt(bckSubtractedCount/bckSubtractedCount**2 + 
                     bckNzDict['NzWeightedSum'][zIndex]/bckNzDict['NzWeightedSum'][zIndex]**2)*delta
    if np.isnan(delta) == True or np.isnan(errDelta) == True:
        print("... delta is nan - i.e., no background galaxies - skipping ...")
        return None
    if errDelta > delta:
        print("... delta highly uncertain (deltaErr > delta) - skipping ...")
        return None
    
    # Optional: de-bias right here (use with caution)
    if zDebias is not None:
        z=z+zDebias*(1+z)
    
    print("... zCluster = %.2f, delta = %.1f, errDelta = %.1f (RADeg = %.6f, decDeg = %.6f) ..." % (z, delta, errDelta, RADeg, decDeg))

    return {'z': z, 'pz': pzWeightedMean, 'zOdds': zOdds, 'pz_z': zArray, 'delta': delta, 'errDelta': errDelta}

#-------------------------------------------------------------------------------------------------------------
def estimateAreaMask(RADeg, decDeg, catalog):
    """Using the objects in the catalog, estimate the area covered - to take care of e.g. survey boundaries
    without needing to get masks for every survey. This works by estimating the average angular separation 
    between sources in the catalog, and using that to draw a circle of that radius around every object. This 
    is then used as the mask in area calculations.
    
    Args:
        RADeg (float): Cluster RA position in decimal degrees.
        decDeg (float): Cluster dec position in decimal degrees.
        catalog (list): List of dictionaries, where each dictionary defines an object in the catalog.
    
    Returns:
        areaMask (2d array): Mask of area covered.
        wcs (astWCS.WCS object): WCS corresponding to areaMask.
    
    """

    # Typical separation in degrees - use to set scale
    RAs=[]
    decs=[]
    for obj in catalog:
        RAs.append(obj['RADeg'])
        decs.append(obj['decDeg'])
    RAs=np.array(RAs)
    decs=np.array(decs)
    catCoords=SkyCoord(ra = RAs, dec = decs, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(catCoords, catCoords, nthneighbor = 2)
    sepDeg=np.median(rDeg.data)
    
    # Make mask - max radius here is from catalogRetriever but should be adjustable elsewhere
    maxRadiusDeg=40/60.0
    widthPix=int(np.ceil((2*maxRadiusDeg)/sepDeg))
    areaMask=np.zeros([widthPix, widthPix])
    im=Image.fromarray(areaMask)
    draw=ImageDraw.Draw(im)
    CRVAL1, CRVAL2=RADeg, decDeg
    xSizeDeg, ySizeDeg=2*maxRadiusDeg, 2*maxRadiusDeg
    xSizePix=float(widthPix)
    ySizePix=float(widthPix)
    xRefPix=xSizePix/2.0
    yRefPix=ySizePix/2.0
    xOutPixScale=xSizeDeg/xSizePix
    yOutPixScale=ySizeDeg/ySizePix
    newHead=pyfits.Header()
    newHead['NAXIS']=2
    newHead['NAXIS1']=xSizePix
    newHead['NAXIS2']=ySizePix
    newHead['CTYPE1']='RA---TAN'
    newHead['CTYPE2']='DEC--TAN'
    newHead['CRVAL1']=CRVAL1
    newHead['CRVAL2']=CRVAL2
    newHead['CRPIX1']=xRefPix+1
    newHead['CRPIX2']=yRefPix+1
    newHead['CDELT1']=-xOutPixScale
    newHead['CDELT2']=xOutPixScale
    newHead['CUNIT1']='DEG'
    newHead['CUNIT2']='DEG'
    wcs=astWCS.WCS(newHead, mode='pyfits')
    for obj in catalog:
        x, y=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
        x=int(x); y=int(y)
        #areaMask[y, x]=1
        draw.ellipse([x-2, y-2, x+2, y+2])
    areaMask=np.array(im)

    return areaMask, wcs
    
#-------------------------------------------------------------------------------------------------------------
def extractArraysFromGalaxyCatalog(catalog, RADeg, decDeg):
    """For convenience: pulls out arrays needed for NGal measurement from galaxy catalog. Needs RADeg, decDeg
    to calculate angular distances from cluster (or random) position.
    
    """
    
    RAs=[]
    decs=[]
    gOdds=[]
    gRedshifts=[]
    #gTemplates=[]
    for obj in catalog:
        RAs.append(obj['RADeg'])
        decs.append(obj['decDeg'])
        if 'pz' in list(obj.keys()):
            gOdds.append(obj['odds'])
            gRedshifts.append(obj['zPhot'])
            #gTemplates.append(obj['bestFitSEDFileName'])
    RAs=np.array(RAs)
    decs=np.array(decs)
    angSepArray=astCoords.calcAngSepDeg(RADeg, decDeg, RAs, decs)
    tanArray=np.radians(angSepArray)
    gOdds=np.array(gOdds)
    gRedshifts=np.array(gRedshifts)
    
    return gOdds, gRedshifts, angSepArray, tanArray
        
#-------------------------------------------------------------------------------------------------------------
def calculateRedshiftAndOdds(pz, zArray, dzOdds = 0.2, method = 'max', zPriorMax = None, zPriorMin = None):
    """Calculates z and BPZ/EAZY style odds for given pz, zArray
    
    Method == 'max' finds z based on max peak in pz
    Method == 'odds' finds z based on z at which odds is maximised
    
    """
    
    # Prior: to avoid worst cases of overestimating z
    # Prior: low-z cut, set by 1 Mpc radius requirement and 9' radius input catalogues (np.radians(9.0/60)*astCalc.da(0.1))
    prior=np.ones(pz.shape)
    if zPriorMax is not None:
        prior[np.greater(zArray, zPriorMax)]=0.0
    if zPriorMin is not None:
        prior[np.less(zArray, zPriorMin)]=0.0
    pz=pz*prior
    norm=np.trapz(pz, zArray)
    pz=pz/norm    
        
    zToIndex_tck=interpolate.splrep(zArray, np.arange(zArray.shape[0]))
    z=zArray[pz.tolist().index(pz.max())]
    
    zOdds=None  # we only calculate this if method == 'odds', but we return always
    
    if method == 'max':
        zMin=z-dzOdds
        zMax=z+dzOdds
        indexMin, indexMax=interpolate.splev([zMin, zMax], zToIndex_tck)
        indexMin=int(round(indexMin))
        indexMax=int(round(indexMax))
        if indexMin < 0:
            indexMin=0
        if indexMax > pz.shape[0]-1:
            indexMax=pz.shape[0]-1
        odds=np.trapz(pz[indexMin:indexMax], zArray[indexMin:indexMax])
    
    elif method == 'odds':
        zOdds=[]
        for z in zArray:
            zMin=z-dzOdds
            zMax=z+dzOdds
            indexMin, indexMax=interpolate.splev([zMin, zMax], zToIndex_tck)
            indexMin=int(round(indexMin))
            indexMax=int(round(indexMax))
            if indexMin < 0:
                indexMin=0
            if indexMax > pz.shape[0]-1:
                indexMax=pz.shape[0]-1
            odds=np.trapz(pz[indexMin:indexMax], zArray[indexMin:indexMax])
            zOdds.append(odds)
        zOdds=np.array(zOdds)
        z=zArray[zOdds.tolist().index(zOdds.max())]
        odds=zOdds[zOdds.tolist().index(zOdds.max())]    
    
    return [z, odds, zOdds]
    
