"""
    Routines for finding cluster redshift from p(z) distributions of galaxies in catalog

    ---
    
    Copyright 2016 Matt Hilton (matt.hilton@mykolab.com)
    
    This file is part of zCluster.

    zCluster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sourcery is distributed in the hope that it will be useful,
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
import pylab as plt
import IPython
plt.matplotlib.interactive(True)

#-------------------------------------------------------------------------------------------------------------
def makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, weightsType, minDistanceMpc = 0.0, maxDistanceMpc = 1.0, 
                   sanityCheckRadiusArcmin = 1.0):
    """Make a N(z) distribution in the direction of the postion (usually of a cluster) given by RADeg, decDeg.
    This is constructed from the p(z) distributions of all galaxies in the catalog (subject to odds cut?),
    with radial weights applied as per weightsType. So, for 'flat1Mpc', the radial weight is 1 if within a
    projected distance of 1 Mpc of RADeg, decDeg, and 0 otherwise.
    
    Since this is aimed at detecting clusters, zPriorMax is used to reign in any spuriously high redshifts
    given our prior knowledge of the depth of a given photometric survey
    
    NOTE: minDistanceMpc only applies to 'flat' weightsType (for local background estimate)
    
    """

    # For speed, we make a big array of all p(z) here
    pzArray=[]
    for obj in catalog:
        if 'pz' in obj.keys():
            pzArray.append(obj['pz'])
    pzArray=np.array(pzArray)
    
    # Make big array of angular separation for each galaxy
    # Get EAZY odds too, in case useful for weighting
    gOdds, gRedshifts, angSepArray, tanArray=extractArraysFromGalaxyCatalog(catalog, RADeg, decDeg)
    
    # Sanity check: if we didn't find a decent number of galaxies close to the cluster centre, 
    # the cluster is probably not actually in the optical catalog footprint (e.g., just outside edge of S82)
    if np.sum(np.less(np.degrees(tanArray), sanityCheckRadiusArcmin/60.0)) == 0:
        print "... no galaxies within %.1f' of cluster position ..." % (sanityCheckRadiusArcmin)
        return None

    # Assume same z grid for all objects, calculate DAs
    for obj in catalog:
        if 'pz_z' in obj.keys():
            zArray=obj['pz_z']
            break
    DAArray=[]
    for z in zArray:
        DAArray.append(astCalc.da(z))
    DAArray=np.array(DAArray)
        
    # NOTE: radial weights should have max value 1, min value 0
    # Koester, Postman et al. style radial weights (projected 2d NFW profile)
    if weightsType == 'NFW':
        rs=0.15 # rs=R200/c c=concentration, presumably - Koester uses 150 kpc
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
        sigmax[mask]=(2*rs*fx[mask])/(x[mask]**2-1)   # ignoring rho_s, which is arbitrary
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
        raise Exception, "didn't understand weightsType"
    
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
    Nz_r=rWeights.sum(axis = 0)                         # Number of galaxies within 1 Mpc at each z (IF all weights were 1)
    Nz_r_odds_total=np.greater(rWeights*gOddsWeights, 0).sum(axis = 0)      # Actual number of galaxies within 1Mpc at each z regardless of weight (if using NFW, radial - for flat, Nz_total == Nz_r)
    Nz_r_odds=np.sum(rWeights*gOddsWeights, axis = 0)   # Number of galaxies within 1 Mpc at each z AND odds > 0.9
    NzWeightedSum=np.sum(pzArray*rWeights*gOddsWeights, axis = 0)
    
    # We may as well keep track of area as well
    area=np.pi*maxDistanceMpc**2 - np.pi*minDistanceMpc**2
    
    # If we want p(z), we can divide pzWeightedSum by Nz_r_odds
    return {'NzWeightedSum': NzWeightedSum, 'Nz': Nz_r_odds, 'Nz_total': Nz_r_odds_total, 'zArray': zArray, 'areaMpc2': area}
    
#-------------------------------------------------------------------------------------------------------------
def estimateClusterRedshift(RADeg, decDeg, catalog, zPriorMin, zPriorMax, weightsType, maxRMpc, 
                            zMethod, sanityCheckRadiusArcmin = 1.0):
    """This does the actual work of estimating cluster photo-z from catalog.
    
    Assumes each object has keys 'pz' (p(z), probability distribution), 'pz_z' (corresponding redshifts at 
    each point in p(z)
            
    """
    
    print "... estimating cluster photo-z ..."

    # Initial sanity check: if we didn't find a decent number of galaxies close to the cluster centre, 
    # the cluster is probably not actually in the optical catalog footprint (e.g., just outside edge of S82)
    gOdds, gRedshifts, angSepArray, tanArray=extractArraysFromGalaxyCatalog(catalog, RADeg, decDeg)
    if np.sum(np.less(np.degrees(tanArray), sanityCheckRadiusArcmin/60.0)) == 0:
        print "... no galaxies within %.1f' of cluster position ..." % (sanityCheckRadiusArcmin)
        return None
        
    # Here, we do a weighted average of the p(z) distribution
    # The weighting is done according to distance from the cluster position at the z corresponding to a
    # given p(z). 
    clusterNzDict=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, weightsType, maxDistanceMpc = 1.0)
    zArray=clusterNzDict['zArray']
    if clusterNzDict['NzWeightedSum'].sum() == 0:
        print "... no galaxies in n(z) - skipping..."
        return None
    
    # This is how we have been calculating redshifts (see estimateClusterRedshift_old below)
    # We only need normFactor to make sure odds is correct
    pzWeightedMean=clusterNzDict['NzWeightedSum']/clusterNzDict['Nz']
    pzWeightedMean[np.isnan(pzWeightedMean)]=0.0
    normFactor=np.trapz(pzWeightedMean, clusterNzDict['zArray'])
    pzWeightedMean=pzWeightedMean/normFactor
    try:
        z, odds, zOdds=calculateRedshiftAndOdds(pzWeightedMean, clusterNzDict['zArray'], dzOdds = 0.05, method = zMethod, zPriorMax = zPriorMax, zPriorMin = zPriorMin)
    except:
        print "Hmm - failed in z, odds calc: what happened?"
        IPython.embed()
        sys.exit()
                
    # Local background - 0.6 deg radius gives us out to 4 Mpc radius at z = 0.1
    zIndex=np.where(zArray == z)[0][0]
    clusterNzDictForSNR=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, 'flat', maxDistanceMpc = maxRMpc)
    localBckNzDict=makeWeightedNz(RADeg, decDeg, catalog, zPriorMax, 'flat', minDistanceMpc = 3.0, maxDistanceMpc = 4.0)
    bckAreaNorm=clusterNzDictForSNR['areaMpc2']/localBckNzDict['areaMpc2']
    bckSubtractedCount=clusterNzDictForSNR['NzWeightedSum'][zIndex]-bckAreaNorm*localBckNzDict['NzWeightedSum'][zIndex]
    delta=bckSubtractedCount/(bckAreaNorm*localBckNzDict['NzWeightedSum'][zIndex])
    errDelta=np.sqrt(bckSubtractedCount/bckSubtractedCount**2 + 
                     localBckNzDict['NzWeightedSum'][zIndex]/localBckNzDict['NzWeightedSum'][zIndex]**2)*delta

    if np.isnan(errDelta) == True:
        errDelta=0
        
    print "... zCluster = %.2f, odds = %.2f, ngal = %d, delta = %.1f, errDelta = %.1f (RADeg = %.6f, decDeg = %.6f) ..." % (z, odds, bckSubtractedCount, delta, errDelta, RADeg, decDeg)

    return {'z': z, 'odds': odds, 'pz': pzWeightedMean, 'zOdds': zOdds, 'pz_z': zArray, 
            'ngal': bckSubtractedCount, 'delta': delta, 'errDelta': errDelta}

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
        if 'pz' in obj.keys():
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
    if zPriorMax != None:
        prior[np.greater(zArray, zPriorMax)]=0.0
    if zPriorMin != None:
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
    
