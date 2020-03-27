"""
    Class for running photo-z code
    
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
import numpy as np
import zCluster
from astLib import *
from scipy import stats
from scipy import interpolate
import string
import glob
import pickle
import time
import IPython

class PhotoRedshiftEngine:
    """A class that calculates galaxy photo-zs, adding photo-z info to a catalog (list of dictionaries)
    in place.
    
    """
    
    def __init__(self, absMagCut, passbandSet = 'SDSS+Ks', zMin = 0.01, zMax = 3.0, zStep = 0.01, 
                 ZPError = 0.0, ZPOffsets = None, templatesDir = None):
        """Sets up the stuff we would otherwise calculate every time, i.e., the templates.
        
        """
                                
        # Redshift grid on which to calculate p(z)
        self.zRange=np.linspace(zMin, zMax, int(((zMax+zStep)-zMin)/zStep))

        # Set up passbands
        passbandsDir=zCluster.__path__[0]+os.path.sep+"passbands"+os.path.sep
        self.passbandsList=[]
        if passbandSet == 'SDSS+Ks':
            self.bands=['u', 'g', 'r', 'i', 'z', 'Ks']
            for band in self.bands:
                if band == 'Ks':
                    self.passbandsList.append(astSED.Passband(passbandsDir+"K_2MASS.res"))
                else:
                    self.passbandsList.append(astSED.Passband(passbandsDir+band+"_SDSS.res"))
        elif passbandSet == 'KiDS-VIKING':
            self.bands=['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']
            for band in self.bands:
                if band in ['Z', 'Y', 'J', 'H', 'Ks']:
                    self.passbandsList.append(astSED.Passband(passbandsDir+"VISTA_Filters_at80K_forETC_%s.dat" % (band), 
                                                              inputUnits = 'nanometres'))
                else:
                    self.passbandsList.append(astSED.Passband(passbandsDir+band+"_SDSS.res"))
        elif passbandSet == 'WIRDS':
            self.bands=['u', 'g', 'r', 'i', 'z', 'J', 'H', 'Ks']
            for band in self.bands:
                if band in ['J', 'H', 'Ks']:
                    self.passbandsList.append(astSED.Passband(passbandsDir+"VISTA_Filters_at80K_forETC_%s.dat" % (band), 
                                                              inputUnits = 'nanometres'))
                else:
                    self.passbandsList.append(astSED.Passband(passbandsDir+band+"_SDSS.res"))
        elif passbandSet == 'PS1':
            self.bands=['g', 'r', 'i', 'z']#, 'y']
            for band in self.bands:
                self.passbandsList.append(astSED.Passband(passbandsDir+band+"PS1.res"))
        elif passbandSet == 'DES':
            self.bands=['g', 'r', 'i', 'z']
            for band in self.bands:
                self.passbandsList.append(astSED.Passband(passbandsDir+band+"DES.res"))
        elif passbandSet == 'DES+WISE':
            # NOTE: WISE passbands from here: http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
            self.bands=['g', 'r', 'z', "w1", "w2"]
            for band in self.bands:
                if band[0] != "w":
                    self.passbandsList.append(astSED.Passband(passbandsDir+band+"DES.res"))
                else:
                    self.passbandsList.append(astSED.Passband(passbandsDir+"RSR-%s.EE.txt" % (band.upper()), inputUnits = "microns"))
        elif passbandSet == 'DECaLS':
            # NOTE: WISE passbands from here: http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
            self.bands=['g', 'r', 'z', "w1", "w2"]
            for band in self.bands:
                if band[0] != "w":
                    self.passbandsList.append(astSED.Passband(passbandsDir+"decam.%s.am1p4.dat.txt" % (band)))
                else:
                    self.passbandsList.append(astSED.Passband(passbandsDir+"RSR-%s.EE.txt" % (band.upper()), inputUnits = "microns"))
        else:
            raise Exception("Unknown passbandSet '%s'" % (passbandSet))
        
        # Allow a global photometric calibration error, in magnitudes, (across all bands)
        # This gets added to photometric errors in catalog in quadrature
        # This helps to get catalogs with extremely small photometric errors to behave
        self.ZPError=ZPError
                
        # This probably doesn't gain us much
        self.effectiveWavelength=[]
        for p in self.passbandsList:
            self.effectiveWavelength.append(p.effectiveWavelength())
        self.effectiveWavelength=np.array(self.effectiveWavelength)
        self.effLMicron2=(self.effectiveWavelength*(1e-10/1e-6))*(self.effectiveWavelength*(1e-10/1e-6))
            
        # SED templates 
        # Mostly lifted from those included with EAZY. CWW is CWW+Kinney, used in BPZ etc..
        # EAZY_v1.0 is a NMF derived (basically a PCA) minimal template set from the PEGASE2.0 SEDs
        # BR07 is a NMF derived (basically a PCA) minimal template set derived from a load of BC03 models, designed
        # to reproduce the colours of galaxies in SDSS.
        self.modelSEDDictList=[]
        if templatesDir is None:
            pickleFileName=None
            self.SEDFiles=glob.glob(zCluster.__path__[0]+os.path.sep+"SED/EAZY_v1.0/*.dat")+glob.glob(zCluster.__path__[0]+os.path.sep+"SED/CWW/*.sed")
            #self.SEDFiles=glob.glob(zCluster.__path__[0]+os.path.sep+"SED/uvista_nmf/*.dat")+glob.glob(zCluster.__path__[0]+os.path.sep+"SED/CWW/*.sed")+glob.glob(zCluster.__path__[0]+os.path.sep+"SED/EAZY_v1.0/*.dat")
            #self.SEDFiles=glob.glob(zCluster.__path__[0]+os.path.sep+"SED/EAZY_v1.1_lines/*.dat")+glob.glob(zCluster.__path__[0]+os.path.sep+"SED/CWW/*.sed")+glob.glob(zCluster.__path__[0]+os.path.sep+"SED/uvista_nmf/*.dat")
            #self.SEDFiles=glob.glob(zCluster.__path__[0]+os.path.sep+"SED/BC03Templates/*.res")#+glob.glob(zCluster.__path__[0]+os.path.sep+"SED/CWW/*.sed")
        else:
            # We'll try a few different extensions
            print(">>> Using custom template set from %s" % (templatesDir))
            pickleFileName=templatesDir+os.path.sep+"templates_%s.pkl" % (passbandSet)
            self.SEDFiles=glob.glob(templatesDir+os.path.sep+"*.res")
            self.SEDFiles=self.SEDFiles+glob.glob(templatesDir+os.path.sep+"*.dat")
            self.SEDFiles=self.SEDFiles+glob.glob(templatesDir+os.path.sep+"*.sed")

        if pickleFileName is not None and os.path.exists(pickleFileName):
            with open(pickleFileName, "rb") as pickleFile:
                unpickler=pickle.Unpickler(pickleFile)
                self.numModels=unpickler.load()
                self.modelFlux=unpickler.load()
                self.templateIndex=unpickler.load()
                self.modelFlux2=self.modelFlux**2
        else:
            self.numModels=len(self.SEDFiles)
            i=0
            t=0
            self.templateIndex=[]
            for f in self.SEDFiles:
                s=astSED.SED()
                s.loadFromFile(f)
                t=t+1
                for z in self.zRange:
                    s.redshift(z)
                    modelSEDDict=s.getSEDDict(self.passbandsList)
                    modelSEDDict['E(B-V)']=None
                    modelSEDDict['ageGyr']=0.0
                    modelSEDDict['z']=z
                    modelSEDDict['fileName']=f 
                    modelSEDDict['modelListIndex']=i
                    modelSEDDict['SED']=s.copy()
                    self.modelSEDDictList.append(modelSEDDict)       
                    self.templateIndex.append(t)
                i=i+1
                del s
            self.templateIndex=np.array(self.templateIndex)                
            modelFlux=[]
            for modelSEDDict in self.modelSEDDictList:
                modelFlux.append(modelSEDDict['flux'])
            self.modelFlux=np.array(modelFlux)  
            self.modelFlux2=self.modelFlux**2
        
        # We can pickle if give a custom templates dir
        if templatesDir is not None:
            with open(pickleFileName, "wb") as pickleFile:
                pickler=pickle.Pickler(pickleFile)
                pickler.dump(self.numModels)
                pickler.dump(self.modelFlux)
                pickler.dump(self.templateIndex)
        
        # These can be fitted for on the fly or just specified as an argument
        self.ZPOffsets=np.zeros(self.modelFlux.shape[1])
        if ZPOffsets is not None:
            self.ZPOffsets=ZPOffsets
            
        # We might use these...
        dlRange=[]
        for z in self.zRange:
            dlRange.append(astCalc.dl(z))
        self.dlRange=np.array(dlRange)
        
        # More sophisticated mag prior... doesn't seem to help
        # Using m* in i-band at z = 0.1 from Popesso et al. and BC03 tau = 0.1 Gyr, solar metallicity, zf = 3
        # Cut is m*-3, i.e., should leave in BCGs
        #MStar=-21.8
        #MStarBand='i'
        #bc03=astSED.BC03Model(zCluster.__path__[0]+os.path.sep+"data"+os.path.sep+"tau0p1Gyr_m62.20")
        #passband=self.passbandsList[self.bands.index(MStarBand)]
        #m=MStar+5.*np.log10(astCalc.dl(0.1)*1e5)
        #magEvo=bc03.getMagEvolution(passband, m, 0.1, 3.0, zStepSize=0.05, magType='AB')   
        #tck=interpolate.splrep(magEvo['z'], magEvo['mag'])
        #self.magPriorCut=interpolate.splev(self.zRange, tck)-5.*np.log10(self.dlRange*1e5)-magsBrighterMStarCut
        #self.magPriorBand=self.bands.index(MStarBand)
        
        # Less sophisticated mag prior...
        self.magPriorCut=absMagCut
        self.magPriorBand=self.bands.index('r')

    def calcZeroPointOffsets(self, galaxyCatalog, zSpecColumn = 'z_spec'):
        """If the galaxyCatalog contains zSpecColumn, fit SEDs at given spec-zs to determine
        zero point offsets. After this is done, they will automatically be applied by calcPhotoRedshifts.
        
        """
        print("... calculating zero point offsets for each band using objects with spec-zs ...")
        ten23=(10**23.0)
        diffMags=[]
        for galaxy in galaxyCatalog:
            if zSpecColumn in galaxy.keys() and galaxy[zSpecColumn] > 0:
                # If we don't have a band, include a ridiculous mag with ridiculous error so zero weight in fit
                magAB=[]
                magErrAB=[]
                for k in self.bands:
                    if k in list(galaxy.keys()):
                        magAB.append(galaxy[k])
                        magErrAB.append(galaxy['%sErr' % (k)])
                    else:
                        magAB.append(99.0)
                        magErrAB.append(99.0)
                magAB=np.array(magAB)
                magErrAB=np.array(magErrAB)
            
                # Incorporate a global photometric calibration error (applied to all bands)
                # Helps to get catalogs with tiny errors to behave
                magErrAB=np.sqrt(magErrAB**2+self.ZPError**2)
                            
                fluxJy=ten23*np.power(10, (-(magAB+48.6)/2.5)) 
                sedFlux=3e-13*fluxJy/self.effLMicron2
                fluxJyErr=ten23*np.power(10, (-(magAB-magErrAB+48.6)/2.5))
                sedFluxErr=(3e-13*fluxJyErr/self.effLMicron2)-sedFlux
                sedFluxErr2=sedFluxErr**2
                # Note we don't actually need to make sedFlux same shape as modelFlux, big speed up, same result
                norm=np.sum((self.modelFlux*sedFlux)/(sedFluxErr2), axis=1)/np.sum(self.modelFlux2/sedFluxErr2, axis=1)
                chiSq=np.sum(((sedFlux-norm.reshape([norm.shape[0], 1])*self.modelFlux)**2)/sedFluxErr2, axis=1)
                chiSq[np.isnan(chiSq)]=1e6   # throw these out, should check this out and handle more gracefully

                chiSqProb=np.exp(-chiSq/2)#stats.chisqprob(chiSq, len(self.bands)-2)
                chiSqProb=chiSqProb.reshape([self.numModels, self.zRange.shape[0]])
                zi=np.argmin(abs(galaxy[zSpecColumn]-self.zRange))
                
                chiSqProb=chiSqProb[:, zi]
                modFlux=(norm.reshape([norm.shape[0], 1])*self.modelFlux)
                modFlux=modFlux.reshape([self.numModels, self.zRange.shape[0], len(magAB)]) 
                modFlux=modFlux[:, zi, :]
                
                modFlux=np.average(modFlux, axis = 0, weights = chiSqProb)
                modFluxJy=(modFlux*self.effLMicron2)/3e-13 
                modMags=-2.5*np.log10(modFluxJy/1e23)-48.6
                diffMags.append(magAB-modMags)
        # Average after taking out differenced 99 values for missing data
        diffMags=np.array(diffMags)
        ZPOffsets=[]
        for i in range(diffMags.shape[1]):
            d=diffMags[:, i]
            ZPOffsets.append(np.median(d[np.less(d, 1)]))
        self.ZPOffsets=np.array(ZPOffsets)
        print("... offsets found: %s mag." % (str(self.ZPOffsets)))

    
    def calcPhotoRedshifts(self, galaxyCatalog, calcMLRedshiftAndOdds = False):
        """Calculates photometric redshifts and adds to the galaxy catalog in place.
        
        NOTE: since normally we're normally only interested in p(z), this only returns
        the maximum likelihood z and BPZ-style odds parameter if calcMLRedshiftAndOdds = True.
        
        """
               
        # Note that we now aren't checking for junk/bad bands in the catalogue
        # That should have been dealt with before getting to this stage
        # And we've also used some more tricks to speed up by factor of 4 or so
        # Things we tried to speed up the below
        # Own trapz doesn't help here (this is probably to do with array size)
        # Bottle neck is chiSqProb calculation (1.1 sec integrated over whole catalog)
        # Equivalent to special.gammaincc((len(self.bands)-2)/2.0, chiSq/2.0)
        # Not found a faster implementation
        ten23=(10**23.0)
        if self.ZPOffsets.sum() != 0:
            print("... applying zero-point offsets:", self.ZPOffsets, "...")
        for galaxy in galaxyCatalog:

            # If we don't have a band, include a ridiculous mag with ridiculous error so zero weight in fit
            magAB=[]
            magErrAB=[]
            for k in self.bands:
                if k in list(galaxy.keys()):
                    magAB.append(galaxy[k])
                    magErrAB.append(galaxy['%sErr' % (k)])
                else:
                    magAB.append(99.0)
                    magErrAB.append(99.0)
            magAB=np.array(magAB)
            magErrAB=np.array(magErrAB)
            
            magAB=magAB-self.ZPOffsets
            
            # Incorporate a global photometric calibration error (applied to all bands)
            # Helps to get catalogs with tiny errors to behave
            magErrAB=np.sqrt(magErrAB**2+self.ZPError**2)
                        
            fluxJy=ten23*np.power(10, (-(magAB+48.6)/2.5)) 
            sedFlux=3e-13*fluxJy/self.effLMicron2
            fluxJyErr=ten23*np.power(10, (-(magAB-magErrAB+48.6)/2.5))
            sedFluxErr=(3e-13*fluxJyErr/self.effLMicron2)-sedFlux
            sedFluxErr2=sedFluxErr**2
            # Note we don't actually need to make sedFlux same shape as modelFlux, big speed up, same result
            norm=np.sum((self.modelFlux*sedFlux)/(sedFluxErr2), axis=1)/np.sum(self.modelFlux2/sedFluxErr2, axis=1)
            chiSq=np.sum(((sedFlux-norm.reshape([norm.shape[0], 1])*self.modelFlux)**2)/sedFluxErr2, axis=1)
            chiSq[np.isnan(chiSq)]=1e6   # throw these out, should check this out and handle more gracefully
            
            #---
            ## Organise model fluxes by z, template number so a[0, 0] == self.modelFlux[0] is template 0 at z = 0
            #a=self.modelFlux.reshape([self.zRange.shape[0], self.numModels, 8]) 
            #print("improve photo-zs: linear combinations?")
            #IPython.embed()
            #sys.exit()
            
            #---
            # This extracts chiSq as function of redshift for the best-fit template only, if we wanted it
            #chiSq=chiSq[self.templateIndex == self.templateIndex[np.argmin(chiSq)]]
            #pz=np.exp(-chiSq/2)
            # This uses all templates at once
            chiSqProb=np.exp(-chiSq/2)#stats.chisqprob(chiSq, len(self.bands)-2)
            chiSqProb=chiSqProb.reshape([self.numModels, self.zRange.shape[0]])
            #pz=np.sum(chiSqProb, axis = 0)
            pz=np.max(chiSqProb, axis = 0)
            # Mag prior
            absMag=magAB[self.magPriorBand]-5.0*np.log10(1e5*self.dlRange)
            pPrior=np.array(np.greater(absMag, self.magPriorCut), dtype = float)
            pz=pz*pPrior
            # Normalise
            pzNorm=np.trapz(pz, self.zRange)
            if pzNorm != 0:
                pz=pz/pzNorm 
            if calcMLRedshiftAndOdds == True:
                zp, odds=self.calculateMLRedshiftAndOdds(pz, dzOdds = 0.2)
                galaxy['zPhot']=zp
                galaxy['odds']=odds
            galaxy['pz']=pz
            galaxy['pz_z']=self.zRange
        t1=time.time()
        

    def calculateMLRedshiftAndOdds(self, pz, dzOdds = 0.2, method = 'max'):
        """Calculates maximum likelihood z and BPZ/EAZY style odds for given pz, zArray
        
        Method == 'max' finds z based on max peak in pz
        Method == 'odds' finds z based on z at which odds is maximised
        
        """
        
        zToIndex_tck=interpolate.splrep(self.zRange, np.arange(self.zRange.shape[0]))
        z=self.zRange[pz.tolist().index(pz.max())]
        
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
            odds=np.trapz(pz[indexMin:indexMax], self.zRange[indexMin:indexMax])
        
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
                odds=np.trapz(pz[indexMin:indexMax], self.zRange[indexMin:indexMax])
                zOdds.append(odds)
            zOdds=np.array(zOdds)
            z=self.zRange[zOdds.tolist().index(zOdds.max())]
            odds=zOdds[zOdds.tolist().index(zOdds.max())]    
    
        return [z, odds]
    
