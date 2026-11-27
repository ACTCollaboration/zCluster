"""

This module contains routines for estimating galaxy cluster photometric redshifts.

"""

import os
import sys
import numpy as np
from astLib import *
from scipy import interpolate
from scipy import stats
from scipy import ndimage
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from PIL import Image
from PIL import ImageDraw
import pylab as plt

#------------------------------------------------------------------------------------------------------------
def getCatalogArrays(catalog):
    """Extracts the arrays needed by the cluster redshift routines from a galaxy catalog (list of dictionaries,
    with photo-zs already calculated). Building these once per galaxy catalog, and passing the result to
    :func:`estimateClusterRedshift`, :func:`makeWeightedNz`, :func:`makeDensityMap` and
    :func:`estimateAreaMask` via their `catalogArrays` argument, saves re-extracting them on every call
    (e.g., in each iteration of the cluster centre finding loop).

    Args:
        catalog (:obj:`list`): Galaxy catalog, as a list of dictionaries.

    Returns:
        A dictionary of arrays, with keys:
            - 'RADeg', 'decDeg': coordinates of all galaxies in the catalog
            - 'pzRADeg', 'pzDecDeg': coordinates of galaxies with p(z) (key 'pz' in the catalog)
            - 'pz': 2d array of p(z) of galaxies with p(z), one row per galaxy
            - 'pz_z': redshift grid corresponding to p(z) (assumed to be the same for all galaxies)
            - 'odds', 'zPhot': odds and maximum likelihood photo-z for galaxies with p(z)

    Note:
        The catalog must not be modified (e.g., by adding or removing galaxies) after this has been
        called, otherwise the arrays will no longer correspond to the catalog.

    """

    hasPZ=['pz' in obj for obj in catalog]
    pzCatalog=[obj for obj, has in zip(catalog, hasPZ) if has]
    arrays={'RADeg': np.array([obj['RADeg'] for obj in catalog]),
            'decDeg': np.array([obj['decDeg'] for obj in catalog]),
            'pzRADeg': np.array([obj['RADeg'] for obj in pzCatalog]),
            'pzDecDeg': np.array([obj['decDeg'] for obj in pzCatalog]),
            'odds': np.array([obj['odds'] for obj in pzCatalog]),
            'zPhot': np.array([obj['zPhot'] for obj in pzCatalog]),
            'pz': np.array([obj['pz'] for obj in pzCatalog]),
            'pz_z': None}
    for obj in catalog:
        if 'pz_z' in obj:
            arrays['pz_z']=obj['pz_z']
            break

    return arrays

#------------------------------------------------------------------------------------------------------------
_DACache={}
def _getDAArray(zArray):
    """Returns angular diameter distance at each redshift in zArray, caching the result (since we normally
    use the same redshift grid every time).

    """
    key=np.asarray(zArray).tobytes()
    if key not in _DACache:
        _DACache[key]=np.array([astCalc.da(z) for z in zArray])
    return _DACache[key]

#------------------------------------------------------------------------------------------------------------
_NFWSplineCache={}
def _getNFWSpline():
    """Returns a spline fit to the projected 2d NFW profile used for 'NFW' radial weights (Koester, Postman
    et al. style), normalised to have max value 1. This is calculated once and then cached.

    """
    if 'tck' not in _NFWSplineCache:
        rs=0.15 # rs=R200/c c=concentration - Koester uses 150 kpc
        r=np.linspace(0.0, 10.0, 2000)
        x=r/rs
        fx=np.zeros(x.shape)
        mask=np.greater(x, 1)
        fx[mask]=1-(2.0/np.sqrt(x[mask]**2-1))*np.arctan(np.sqrt((x[mask]-1)/(x[mask]+1)))
        mask=np.logical_and(np.greater(x, 0), np.less(x, 1))
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
        _NFWSplineCache['tck']=interpolate.splrep(r, sigmax)
    return _NFWSplineCache['tck']

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

#---------------------------------------------------------------------------------------------------
def makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, maxDistDegrees):
    """Fills (in place) the 2d array degreesMap with distance in degrees from the given position,
    out to some user-specified maximum distance.

    Args:
        degreesMap (:obj:`np.ndarray`): Map (2d array) that will be filled with angular distance
            from the given coordinates. Probably you should feed in an array set to some extreme
            initial value (e.g., 1e6 everywhere) to make it easy to filter for pixels near the
            object coords afterwards.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to degreesMap.
        RADeg (float): RA in decimal degrees of position of interest (e.g., object location).
        decDeg (float): Declination in decimal degrees of position of interest (e.g., object
            location).
        maxDistDegrees: The maximum radius out to which distance will be calculated.

    Returns:
        A map (2d array) of distance in degrees from the given position,
        (min x, max x) pixel coords corresponding to maxDistDegrees box,
        (min y, max y) pixel coords corresponding to maxDistDegrees box

    Note:
        This routine measures the pixel scale local to the given position, then assumes that it
        does not change. So, this routine may only be accurate close to the given position,
        depending upon the WCS projection used.

    """

    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)

    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=degreesMap.shape[0]
    X=degreesMap.shape[1]

    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX > X:
        maxX=X
    if minY < 0:
        minY=0
    if maxY > Y:
        maxY=Y

    xDeg=(np.arange(degreesMap.shape[1])-x0)*xPixScale
    yDeg=(np.arange(degreesMap.shape[0])-y0)*yPixScale
    for i in range(minY, maxY):
        degreesMap[i][minX:maxX]=np.sqrt(yDeg[i]**2+xDeg[minX:maxX]**2)

    return degreesMap, [minX, maxX], [minY, maxY]

#------------------------------------------------------------------------------------------------------------
def makeBlankMap(RADeg, decDeg, sizePix, sizeDeg):
    """Makes a square blank map with a WCS.
    
    Returns:
        - mapData (2d array)
        - wcs (astWCS.WCS object)
    
    """
    CRVAL1, CRVAL2=RADeg, decDeg
    xSizeDeg, ySizeDeg=sizeDeg, sizeDeg
    xSizePix=float(sizePix)
    ySizePix=float(sizePix)
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
    newHead['CDELT2']=xOutPixScale    # Makes more sense to use same pix scale
    newHead['CUNIT1']='deg'
    newHead['CUNIT2']='deg'
    wcs=astWCS.WCS(newHead, mode='pyfits')
    
    return np.zeros([int(sizePix), int(sizePix)], dtype = float), wcs

#-------------------------------------------------------------------------------------------------------------
def makeDensityMap(RADeg, decDeg, catalog, z, dz = 0.1, rMaxMpc = 1.5, sizeMpc = 8.0, MpcPerPix = 0.1,
                   gaussSmoothPix = 2, outFITSFileName = None, outPlotFileName = None, catalogArrays = None):
    """Makes a projected density map within +/- dz of the given redshift, smoothed using a Gaussian kernel,
    then calculates the centroid location, centroid shift (with respect to the original position), and an
    asymmetry statistic.
    
    Args:
        RADeg (float): RA coordinate (decimal degrees) of the center of the output map.
        decDeg (float): dec coordinate (decimal degrees) of the center of the output map.
        z (float): Central redshift at which the projected density will be calculated.
        dz (float, optional): Sets the integration range z +/- dz over which the projected density is
            calculated.
        rMaxMpc (float, optional): The maximum radial distance in Mpc (from RADeg, decDeg) in which the
            centroid search is performed.
        sizeMpc (float, optional): The size of the density map (which is square) in projected Mpc.
        MpcPerPix (float, optional): The pixel size, in projected Mpc, for the output density map.
        gaussSmoothPix (int, optional): The size of the Gaussian smoothing kernel applied to the output
            projected density map, in pixels.
        outFITSFileName (str, optional): If given, write the density map as a FITS file to the given
            location.
        makePlot (str, optional): If given, write a plot of the density map to the given location (the
            file format is determined by the extension, and supported formats depend on the matplotlib
            backend used).
        catalogArrays (dict, optional): Output of :func:`getCatalogArrays` for this catalog. If not given,
            this will be calculated.

    Returns:
        A dictionary with keys:
            - 'map': the output projected density map (2d array)
            - 'wcs': the output WCS
            - 'cRADeg': the RA of the centroid (peak) in the density map
            - 'cDecDeg': the dec of the centroid (peak) in the density map
            - 'offsetArcmin': size of the offset of the centroid (arcmin) from the original RA, dec
                coords
            - 'offsetMpc': size of the offset of the centroid (projected Mpc) from the original RA,
                dec coordinates
            - 'CS': centroid shift parameter (Mpc)
            - 'A': asymmetry parameter
    
    """
    
    # We go out to 4 Mpc for delta
    DA=astCalc.da(z)
    sizeDeg=np.degrees(sizeMpc/DA)
    sizePix=int(sizeMpc/MpcPerPix)
    d, wcs=makeBlankMap(RADeg, decDeg, sizePix, sizeDeg)
    if catalogArrays is None:
        catalogArrays=getCatalogArrays(catalog)
    if catalogArrays['pz'].shape[0] > 0:
        pixCoords=np.array(wcs.wcs2pix(catalogArrays['pzRADeg'], catalogArrays['pzDecDeg'])).reshape(-1, 2)
        x=np.round(pixCoords[:, 0]).astype(int)
        y=np.round(pixCoords[:, 1]).astype(int)
        zMask=np.logical_and(catalogArrays['pz_z'] > z-dz, catalogArrays['pz_z'] < z+dz)
        v=np.trapezoid(catalogArrays['pz'][:, zMask], catalogArrays['pz_z'][zMask], axis = 1)
        inMap=np.logical_and(np.logical_and(x >= 0, x < d.shape[1]), np.logical_and(y >= 0, y < d.shape[0]))
        np.add.at(d, (y[inMap], x[inMap]), v[inMap])
    d=ndimage.gaussian_filter(d, gaussSmoothPix)
    
    # Centre shift
    xpos,ypos=np.where(d==d.max())
    CS=np.sqrt((xpos-d.shape[1]/2)**2+(ypos-d.shape[0]/2)**2)[0]
    CSinMpc=CS*MpcPerPix

    # Asymmetry parameter
    d1=np.zeros(d.shape)+d
    d1=np.flipud(np.fliplr(d1))
    Adiff=np.power(d-d1, 2)
    A=np.sqrt((np.power(d-d1, 2).sum())/(2*np.power(d, 2).sum()))
    
    # We restrict centre finding to within 1.5 Mpc radius
    ys=np.array([np.arange(d.shape[0])]*d.shape[1]).transpose()
    xs=np.array([np.arange(d.shape[1])]*d.shape[0])
    coords=np.array(wcs.pix2wcs(xs.flatten(), ys.flatten()))
    ras=coords[:, 0]
    decs=coords[:, 1]
    rDegMap=astCoords.calcAngSepDeg(RADeg, decDeg, ras, decs).reshape(d.shape)
    rMpcMap=np.radians(rDegMap)*DA
    rMask=np.less(rMpcMap, rMaxMpc)
    y, x=np.where((d*rMask) == (d*rMask).max())
    cRADeg, cDecDeg=wcs.pix2wcs(x[0], y[0])
    offsetArcmin=astCoords.calcAngSepDeg(cRADeg, cDecDeg, RADeg, decDeg)*60
    offsetMpc=np.radians(offsetArcmin/60.)*DA
    
    if outFITSFileName is not None:
        astImages.saveFITS(outFITSFileName, d, wcs)

    if outPlotFileName is not None:
        axes=[0, 0, 1, 1]
        axesLabels="sexagesimal"    # Avoid dealing with axis flips
        figSize=(8.25, 8.25)
        fig=plt.figure(figsize = figSize, dpi = 300)
        p=astPlots.ImagePlot(d, wcs, cutLevels = [d.min(), d.max()],
                             colorMapName = 'magma', axes = axes,
                             axesLabels = axesLabels)

        if sizeMpc > 40:
            scaleBarMpc=5
        elif sizeMpc > 10 and sizeMpc < 40:
            scaleBarMpc=2
        elif sizeMpc < 10:
            scaleBarMpc=1
        scaleBarMpc=5#int(round(sizeMpc/5))
        scaleBarSizeArcmin=np.degrees(scaleBarMpc/astCalc.da(z))*60
        p.addScaleBar('SE', scaleBarSizeArcmin*60.0, color='cyan', fontSize=16, width=2.0, label = "%d Mpc" % (scaleBarMpc))
        plt.savefig(outPlotFileName)

    return {'map': d, 'wcs': wcs, 'cRADeg': cRADeg, 'cDecDeg': cDecDeg, 
            'offsetArcmin': offsetArcmin, 'offsetMpc': offsetMpc, 'CS': CS, 'A': A}

#-------------------------------------------------------------------------------------------------------------
_areaMapCache={}
def _getAreaMapsForPosition(areaMask, wcs, RADeg, decDeg):
    """Returns (cropped) maps of distance in degrees from RADeg, decDeg and of the area in square degrees of
    each pixel in areaMask, as used by :func:`makeWeightedNz`. The most recent result is cached (keyed on the
    identity of areaMask and wcs, and the position), as makeWeightedNz is called more than once at the same
    position when estimating a cluster redshift. The pixel area map does not depend on position, so it is
    cached separately.

    """
    cached=_areaMapCache.get('pixArea')
    if cached is None or cached[0] is not areaMask or cached[1] is not wcs:
        _areaMapCache['pixArea']=(areaMask, wcs, getPixelAreaDeg2Map(areaMask, wcs))
        _areaMapCache.pop('position', None)
    pixAreaMap=_areaMapCache['pixArea'][2]
    cached=_areaMapCache.get('position')
    if cached is None or cached[0] != (RADeg, decDeg):
        rDegMap=np.zeros(areaMask.shape)
        rDegMap, xRange, yRange=makeDegreesDistanceMap(rDegMap, wcs, RADeg, decDeg, 2.0)
        rDegMap=rDegMap[yRange[0]:yRange[1], xRange[0]:xRange[1]]
        areaMap=areaMask[yRange[0]:yRange[1], xRange[0]:xRange[1]]*pixAreaMap[yRange[0]:yRange[1], xRange[0]:xRange[1]]
        _areaMapCache['position']=((RADeg, decDeg), rDegMap, areaMap)
    return _areaMapCache['position'][1], _areaMapCache['position'][2]

#-------------------------------------------------------------------------------------------------------------
def makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, weightsType, minDistanceMpc = 0.0, maxDistanceMpc = 1.0, 
                   applySanityCheckRadius = True, sanityCheckRadiusArcmin = 1.0, areaMask = None, wcs = None,
                   catalogArrays = None):
    """Make a N(z) distribution in the direction of the postion (usually of a cluster) given by RADeg, decDeg.
    This is constructed from the p(z) distributions of all galaxies in the catalog, with radial weights applied
    as per weightsType. So, for 'flat1Mpc', the radial weight is 1 if within a projected distance of 1 Mpc of 
    RADeg, decDeg, and 0 otherwise.
    
    zPriorMax can be used to reign in any spuriously high redshifts given our prior knowledge of the depth of a
    given photometric survey. However, this is now largely redundant as the mag-based prior does most of the
    work (see PhotoRedshiftEngine.py).
    
    catalogArrays is the output of :func:`getCatalogArrays` for the catalog - if not given, it will be
    calculated.

    NOTE: minDistanceMpc only applies to 'flat' weightsType (for local background estimate)
    
    """

    if catalogArrays is None:
        catalogArrays=getCatalogArrays(catalog)
    pzArray=catalogArrays['pz']
    gOdds=catalogArrays['odds']

    # Angular separation of each galaxy from the given position
    angSepArray=astCoords.calcAngSepDeg(RADeg, decDeg, catalogArrays['RADeg'], catalogArrays['decDeg'])
    tanArray=np.radians(angSepArray)
    
    # Sanity check: if we didn't find a decent number of galaxies close to the cluster centre, 
    # the cluster is probably not actually in the optical catalog footprint (e.g., just outside edge of S82)
    if applySanityCheckRadius == True and np.sum(np.less(np.degrees(tanArray), sanityCheckRadiusArcmin/60.0)) == 0:
        print("... no galaxies within %.1f' of cluster position ..." % (sanityCheckRadiusArcmin))
        return None

    # Assume same z grid for all objects, calculate DAs
    zArray=catalogArrays['pz_z']
    DAArray=_getDAArray(zArray)

    # Simplest odds weighting scheme is to cut everything with odds < some value
    # Whatever happens, we definitely want to cut anything with odds == 0, as its z will default to zMin
    # NOTE: leave set at gOddsCut=0.5 - this fixed null test problems and didn't negatively impact z accuracy
    # Galaxies with odds < gOddsCut have zero weight in everything below, so we just leave them out
    gOddsCut=0.5
    keep=np.logical_not(np.less(gOdds, gOddsCut))
    #keep[:]=True # uncomment to disable gOdds weights
    pzArray=pzArray[keep]
    rMpc=np.outer(tanArray[keep], DAArray)
        
    # NOTE: radial weights should have max value 1, min value 0
    if weightsType == 'NFW':
        # Koester, Postman et al. style radial weights (projected 2d NFW profile)
        # We have to match area to background: so, zero weight all galaxies outside maxDistanceMpc radius
        rWeights=np.zeros(rMpc.shape)
        inRadius=np.logical_not(np.greater(rMpc, maxDistanceMpc))
        rWeights[inRadius]=interpolate.splev(rMpc[inRadius], _getNFWSpline())
        rWeights[np.less(rWeights, 0)]=0.0
    
    elif weightsType == 'flat':
        # Calculate weights - set to zero for objects outside of fiducial projected distance
        # 1 inside 1 Mpc radius, 0 otherwise
        # Use this for background (the 1 Mpc radius just normalises area
        rWeights=np.ones(rMpc.shape)
        rWeights[np.greater(rMpc, maxDistanceMpc)]=0.0
        rWeights[np.less(rMpc, minDistanceMpc)]=0.0
        
    elif weightsType == 'radial':
        # Calculate weights - set to zero for objects outside of fiducial projected distance
        # 1 if at centre coords, falls linearly with increasing radius
        rWeights=1.0-rMpc/maxDistanceMpc
        rWeights[np.greater(rMpc, maxDistanceMpc)]=0.0
    
    else:
        raise Exception("didn't understand weightsType")
    
    # NOTE: zPrior applied in calculating cluster redshift later, not in counting galaxies here
    Nz_r_odds_total=np.greater(rWeights, 0).sum(axis = 0)   # Actual number of galaxies within 1Mpc at each z AND odds > cut, regardless of weight (if using NFW, radial - for flat, Nz_total == Nz_r)
    Nz_r_odds=np.sum(rWeights, axis = 0)                    # Number of galaxies within 1 Mpc at each z AND odds > cut
    NzWeightedSum=np.sum(pzArray*rWeights, axis = 0)
    
    # We may as well keep track of area as well
    if areaMask is None:
        areaMpc2=np.array([np.pi*maxDistanceMpc**2 - np.pi*minDistanceMpc**2]*len(zArray))
    else:
        rDegMap, areaMap=_getAreaMapsForPosition(areaMask, wcs, RADeg, decDeg)
        areaMpc2=[]
        for DA in DAArray:
            areaScaling=1/np.power(np.degrees(np.arctan(1.0/DA)), 2)
            max_rDegCut=np.degrees(np.arctan(maxDistanceMpc/DA))
            min_rDegCut=np.degrees(np.arctan(minDistanceMpc/DA))
            areaDeg2=areaMap[np.logical_and(rDegMap >= min_rDegCut, rDegMap < max_rDegCut)].sum()
            areaMpc2.append(areaScaling*areaDeg2)
        areaMpc2=np.array(areaMpc2)

    # If we want p(z), we can divide pzWeightedSum by Nz_r_odds
    return {'NzWeightedSum': NzWeightedSum, 'Nz': Nz_r_odds, 'Nz_total': Nz_r_odds_total, 'zArray': zArray, 
            'DAArray': DAArray, 'areaMpc2': areaMpc2}
    
#-------------------------------------------------------------------------------------------------------------
def estimateClusterRedshift(RADeg, decDeg, catalog, zPriorMin, zPriorMax, weightsType, maxRMpc, 
                            zMethod, maskMap = None, maskWCS = None, sanityCheckRadiusArcmin = 1.0, 
                            bckCatalog = [], bckAreaDeg2 = None, filterDeltaValues = True,
                            minBackgroundAreaMpc2 = 11.0, catalogArrays = None, bckCatalogArrays = None):
    """This does the actual work of estimating cluster photo-z from catalog.
    
    Assumes each object has keys 'pz' (p(z), probability distribution), 'pz_z' (corresponding redshifts at 
    each point in p(z)
    
    If bckCatalog and bckAreaDeg2 are given, then these are used for the background estimation. If not, then
    a projected 3-4 Mpc annulus is used.
    
    If zDebias is not None, the final output redshift ('z') will be:
        
        z_final = z_initial + zDebias*(1+z_initial) 
    
    This is obviously a fudge, so use with caution (only used for some surveys - see bin/zCluster).

    catalogArrays and bckCatalogArrays are the output of :func:`getCatalogArrays` for catalog and bckCatalog
    respectively. If not given, they will be calculated. When calling this routine more than once for the
    same catalog (e.g., when iterating on the cluster position), it is quicker to calculate these once and
    pass them in.
            
    """

    if catalogArrays is None:
        catalogArrays=getCatalogArrays(catalog)
    if len(bckCatalog) > 0 and bckCatalogArrays is None:
        bckCatalogArrays=getCatalogArrays(bckCatalog)

    # Initial sanity check: if we didn't find a decent number of galaxies close to the cluster centre, 
    # the cluster is probably not actually in the optical catalog footprint (e.g., just outside edge of S82)
    angSepArray=astCoords.calcAngSepDeg(RADeg, decDeg, catalogArrays['RADeg'], catalogArrays['decDeg'])
    tanArray=np.radians(angSepArray)
    if np.sum(np.less(np.degrees(tanArray), sanityCheckRadiusArcmin/60.0)) == 0:
        return None

    # Automated area mask calculation, using only catalogs
    # Probably not super accurate but better than nothing - will at least take care of survey boundaries
    if maskMap is None and maskWCS is None:
        areaMask, wcs=estimateAreaMask(RADeg, decDeg, catalog, catalogArrays = catalogArrays)
    else:
        areaMask, wcs=maskMap, maskWCS

    # Here, we do a weighted average of the p(z) distribution
    # The weighting is done according to distance from the cluster position at the z corresponding to a
    # given p(z). 
    clusterNzDict=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, weightsType, maxDistanceMpc = 1.0,
                                 areaMask = areaMask, wcs = wcs, catalogArrays = catalogArrays)
    zArray=clusterNzDict['zArray']
    if clusterNzDict['NzWeightedSum'].sum() == 0:
        return None

    # This is how we have been calculating redshifts
    # We only need normFactor to make sure odds is correct
    pzWeightedMean=clusterNzDict['NzWeightedSum']/clusterNzDict['Nz']
    pzWeightedMean[np.isnan(pzWeightedMean)]=0.0
    normFactor=np.trapezoid(pzWeightedMean, clusterNzDict['zArray'])
    pzWeightedMean=pzWeightedMean/normFactor
    z, odds, zOdds=calculateRedshiftAndOdds(pzWeightedMean, clusterNzDict['zArray'], dzOdds = 0.05, method = zMethod,
                                            zPriorMax = zPriorMax, zPriorMin = zPriorMin)

    # Background used in delta (density contrast) measurement
    # We also re-do the N(z) measurement in the cluster region so that this makes sense
    clusterNzDictForSNR=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, 'flat', maxDistanceMpc = maxRMpc,
                                       catalogArrays = catalogArrays)
    if len(bckCatalog) > 0:
        # Global
        bckNzDict=makeWeightedNz(RADeg, decDeg, bckCatalog, zPriorMax, 'flat', minDistanceMpc = 3.0, maxDistanceMpc = 1e9, 
                                 applySanityCheckRadius = False, catalogArrays = bckCatalogArrays)
        bckNzDict['areaMpc2']=((np.radians(np.sqrt(bckAreaDeg2))*bckNzDict['DAArray'])**2)
    else:
        # Local background - 0.6 deg radius gives us out to 4 Mpc radius at z = 0.1
        bckNzDict=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, 'flat', minDistanceMpc = 3.0, maxDistanceMpc = 4.0,
                                 areaMask = areaMask, wcs = wcs, catalogArrays = catalogArrays)
    # Masking sanity check that isn't fully implemented (in some cases would just need to bump up mask resolution)
    if max(bckNzDict['areaMpc2']) > (np.pi*(4**2-3**2)): 
        print("... WARNING: max(areaMpc2) = %.3f > %.3f at (RADeg, decDeg) = (%.6f, %.6f)" % (max(bckNzDict['areaMpc2']), np.pi*(4**2-3**2), RADeg, decDeg))        
        #raise Exception("Mask resolution is too coarse - area from mask > area of circular annulus") 
    
    # Delta calculation with bootstrap error over the whole z range
    # For the first part of the mask, we allow background area == 1/2 of the 3-4 Mpc ring (previous min value was 20)
    # This allows us to get out to survey boundaries but will have large delta error bars
    validMask=np.greater_equal(bckNzDict['areaMpc2'], minBackgroundAreaMpc2) 
    bckAreaNorm=np.zeros(len(clusterNzDictForSNR['areaMpc2']))
    bckAreaNorm[validMask]=clusterNzDictForSNR['areaMpc2'][validMask]/bckNzDict['areaMpc2'][validMask]
    validMask=np.logical_and(validMask, clusterNzDictForSNR['NzWeightedSum'] > 0)
    validMask=np.logical_and(validMask, bckNzDict['NzWeightedSum'] > 0)
    if validMask.sum() > 0:
        nc=clusterNzDictForSNR['NzWeightedSum'][validMask]
        nb=bckNzDict['NzWeightedSum'][validMask]
        err_nc=np.sqrt(nc)
        err_nb=np.sqrt(nb)
        A=bckAreaNorm[validMask]
        delta=(nc/(A*nb))-1
        # 5000 bootstrap samples needed for this to converge at 2 decimal places after rounding
        # NOTE: Drawn as one array, with each row being [nc draws, nb draws] - this consumes the random number
        # stream in the same order as drawing poisson(nc), poisson(nb) alternately in a loop
        numBootstraps=5000
        lam=np.broadcast_to(np.concatenate([nc, nb]), (numBootstraps, 2*len(nc)))
        bs=np.random.poisson(lam)
        bs_nc=bs[:, :len(nc)]
        bs_nb=bs[:, len(nc):]
        bs_nb[bs_nb == 0]=1 # Just to avoid div 0 warnings... background should never be 0 anyway
        bs_delta=(bs_nc/(A*bs_nb))-1
        errDelta=np.std(bs_delta, axis = 0)
        errDelta=np.round(errDelta, 2)
        
        # This uses the old delta method (effectively)
        zDelta=zArray[validMask]
        try:
            deltaIndex=np.where(zDelta == z)[0][0]
            delta_at_z=delta[deltaIndex]
            errDelta_at_z=errDelta[deltaIndex]
        except:
            print("... no valid delta value ...")
            return None
        
        # Optionally filter zOdds according to whether delta is > 3 sigma
        if filterDeltaValues == True:
            if zMethod == 'odds':
                pzToCheck=zOdds
            elif zMethod == 'max':
                pzToCheck=pzWeightedMean
            pzToCheck=applyUniformPrior(pzToCheck, zArray, zPriorMax = zPriorMax, zPriorMin = zPriorMin) 
            pzToCheck_validDelta=np.greater(delta/errDelta, 3)*pzToCheck[validMask]
            zIndex=np.argmax(pzToCheck_validDelta)
            z=zDelta[zIndex]
            delta_at_z=delta[zIndex]
            errDelta_at_z=errDelta[zIndex]
            if z < zPriorMin or z > zPriorMax:
                print("... not within z prior range ...")
                return None
        assert(np.isinf(delta_at_z) == False)
            
    else:
        print("... background area too small - skipping ...")
        return None
    
    return {'z': z, 'pz': pzWeightedMean, 'zOdds': zOdds, 'pz_z': zArray, 'delta': delta_at_z, 'errDelta': errDelta_at_z,
            'areaMask': areaMask, 'wcs': wcs}

#-------------------------------------------------------------------------------------------------------------
def estimateAreaMask(RADeg, decDeg, catalog, catalogArrays = None):
    """Using the objects in the catalog, estimate the area covered - to take care of e.g. survey boundaries
    without needing to get masks for every survey. This works by estimating the average angular separation 
    between sources in the catalog, and using that to draw a circle of that radius around every object. This 
    is then used as the mask in area calculations.
    
    Args:
        RADeg (float): Cluster RA position in decimal degrees.
        decDeg (float): Cluster dec position in decimal degrees.
        catalog (list): List of dictionaries, where each dictionary defines an object in the catalog.
        catalogArrays (dict, optional): Output of :func:`getCatalogArrays` for this catalog. If not given,
            this will be calculated.
    
    Returns:
        areaMask (2d array): Mask of area covered.
        wcs (astWCS.WCS object): WCS corresponding to areaMask.
    
    """

    if catalogArrays is None:
        catalogArrays=getCatalogArrays(catalog)
    RAs=catalogArrays['RADeg']
    decs=catalogArrays['decDeg']

    # Typical separation in degrees - use to set scale
    # This doesn't depend on position, so we cache it in catalogArrays
    if 'medianSepDeg' not in catalogArrays:
        catCoords=SkyCoord(ra = RAs, dec = decs, unit = 'deg')
        xIndices, rDeg, sep3d = match_coordinates_sky(catCoords, catCoords, nthneighbor = 2)
        catalogArrays['medianSepDeg']=np.median(rDeg.data)
    sepDeg=catalogArrays['medianSepDeg']
    
    # Make mask - max radius here is from catalogRetriever but should be adjustable elsewhere
    maxRadiusDeg=40/60.0
    widthPix=int(np.ceil((2*maxRadiusDeg)/sepDeg))
    areaMask=np.zeros([widthPix, widthPix])
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
    newHead['CUNIT1']='deg'
    newHead['CUNIT2']='deg'
    wcs=astWCS.WCS(newHead, mode='pyfits')
    # Mark a small circle around every object. Every circle is the same shape, so we draw one once (with PIL,
    # as we used to for every object), and then stamp its pixels onto the mask at each object position
    stampHalfWidth=2
    stampIm=Image.fromarray(np.zeros([2*stampHalfWidth+1, 2*stampHalfWidth+1]))
    ImageDraw.Draw(stampIm).ellipse([0, 0, 2*stampHalfWidth, 2*stampHalfWidth])
    stamp=np.array(stampIm)
    dy, dx=np.nonzero(stamp)
    inkValue=stamp[dy[0], dx[0]]
    if len(RAs) > 0:
        pixCoords=np.array(wcs.wcs2pix(RAs, decs)).reshape(-1, 2)
        x=pixCoords[:, 0].astype(int)
        y=pixCoords[:, 1].astype(int)
        ys=(y.reshape(-1, 1)-stampHalfWidth+dy).ravel()
        xs=(x.reshape(-1, 1)-stampHalfWidth+dx).ravel()
        inMask=np.logical_and(np.logical_and(xs >= 0, xs < widthPix), np.logical_and(ys >= 0, ys < widthPix))
        areaMask[ys[inMask], xs[inMask]]=inkValue
    areaMask=areaMask.astype(stamp.dtype)

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
def applyUniformPrior(pz, zArray, zPriorMin = None, zPriorMax = None):
    """Applies uniform redshift prior to pz between zPriorMin and zPriorMax.
    
    """
    
    prior=np.ones(pz.shape)
    if zPriorMax is not None:
        prior[np.greater(zArray, zPriorMax)]=0.0
    if zPriorMin is not None:
        prior[np.less(zArray, zPriorMin)]=0.0
    pz=pz*prior
    norm=np.trapezoid(pz, zArray)
    if norm != 0:
        pz=pz/norm
    else:
        pz[:]=0
    
    return pz
    
#-------------------------------------------------------------------------------------------------------------
def calculateRedshiftAndOdds(pz, zArray, dzOdds = 0.2, method = 'max', zPriorMax = None, zPriorMin = None):
    """Calculates z and BPZ/EAZY style odds for given pz, zArray
    
    Method == 'max' finds z based on max peak in pz
    Method == 'odds' finds z based on z at which odds is maximised
    
    """
        
    # Prior: to avoid worst cases of overestimating z
    pz=applyUniformPrior(pz, zArray, zPriorMax = zPriorMax, zPriorMin = zPriorMin) 
        
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
        odds=np.trapezoid(pz[indexMin:indexMax], zArray[indexMin:indexMax])
        
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
            odds=np.trapezoid(pz[indexMin:indexMax], zArray[indexMin:indexMax])
            zOdds.append(odds)
        zOdds=np.array(zOdds)
        z=zArray[np.argmax(zOdds)]
        odds=zOdds[np.argmax(zOdds)]
        #--
        norm=np.trapezoid(zOdds, zArray)
        zOdds=zOdds/norm
            
    return [z, odds, zOdds]
    
