# Cython routines taken from nemo

#include "python.pxi"
#include "numpy.pxi"

from astLib import *
import numpy as np
cimport numpy as np
import cython
import math
import time
import sys
from cython cimport parallel

#-------------------------------------------------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def makeDegreesDistanceMap(np.ndarray[np.float64_t, ndim=2] degreesMap, wcs, RADeg, decDeg, maxDistDegrees):
    """Fills (in place) the 2d array degreesMap with distance in degrees from the given position, out to some
    user-specified maximum distance.
    
    Args:
        degreesMap (:obj:`np.ndarray`): Map (2d array) that will be filled with angular distance from the 
            given coordinates. Probably you should feed in an array set to some extreme initial value (e.g.,
            1e6 everywhere) to make it easy to filter for pixels near the object coords afterwards.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to degreesMap.
        RADeg (float): RA in decimal degrees of position of interest (e.g., object location).
        decDeg (float): Declination in decimal degrees of position of interest (e.g., object location).
        maxDistDegrees: The maximum radius out to which distance will be calculated.
    
    Returns:
        A map (2d array) of distance in degrees from the given position,
        (min x, max x) pixel coords corresponding to maxDistDegrees box, 
        (min y, max y) pixel coords corresponding to maxDistDegrees box
    
    Note:
        This routine measures the pixel scale local to the given position, then assumes that it does not 
        change. So, this routine may only be accurate close to the given position, depending upon the WCS
        projection used.
    
    """
    
    # Pixel distance grid            
    cdef float x0, y0, ra0, dec0, ra1, dec1, xPixScale, yPixScale
    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY, xDistPix, yDistPix
    #cdef np.ndarray[np.float64_t, ndim=2] degreesMap
    
    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    
    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=degreesMap.shape[0]
    X=degreesMap.shape[1]
        
    # Real space map of angular distance in degrees, but only consider values near x0, y0
    #degreesMap=np.ones([Y, X], dtype=np.float64)*1e6 # Allocating this was the bottleneck
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
    for y in range(minY, maxY):
        for x in range(minX, maxX):
            yDeg=(y-y0)*yPixScale
            xDeg=(x-x0)*xPixScale
            degreesMap[y, x]=math.sqrt(xDeg*xDeg+yDeg*yDeg)

    return degreesMap, [minX, maxX], [minY, maxY]

#-------------------------------------------------------------------------------------------------------------
def makeXYDegreesDistanceMaps(np.ndarray[np.float64_t, ndim=2] data, wcs, RADeg, decDeg, maxDistDegrees):
    """Returns an array of distance along x, y axes in degrees from given position. maxDistDegrees sets the 
    limit around the given RADeg, decDeg position to which the distance is calculated.
    
    """
    
    # Pixel distance grid            
    cdef float x0, y0, ra0, dec0, ra1, dec1, xPixScale, yPixScale
    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY, xDistPix, yDistPix
    cdef np.ndarray[np.float64_t, ndim=2] xDegreesMap
    cdef np.ndarray[np.float64_t, ndim=2] yDegreesMap

    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    
    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=data.shape[0]
    X=data.shape[1]
    
    # Real space map of angular distance in degrees, but only consider values near x0, y0
    xDegreesMap=np.ones([Y, X], dtype=np.float64)*1e6
    yDegreesMap=np.ones([Y, X], dtype=np.float64)*1e6
    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX >= X:
        maxX=X-1
    if minY < 0:
        minY=0
    if maxY >= Y:
        maxY=Y-1
    for y in range(minY, maxY):
        for x in range(minX, maxX):
            yDeg=(y-y0)*yPixScale
            xDeg=(x-x0)*xPixScale
            xDegreesMap[y, x]=xDeg
            yDegreesMap[y, x]=yDeg
    
    return [xDegreesMap, yDegreesMap]
