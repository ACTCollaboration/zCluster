"""

This module contains a class for estimating galaxy photometric redshifts using a template fitting method.

"""

import os
import numpy as np
import zCluster
from zCluster import stellarmass as sm
from astLib import *
from scipy import stats
from scipy import interpolate
import string
import glob
import pickle
import time

class PhotoRedshiftEngine:
    """A class that calculates galaxy photo-zs, adding photo-z info to a catalog (list of dictionaries)
    in place.
    
    """
    
    def __init__(self, absMagCut, passbandSet = 'SDSS+Ks', zMin = 0.01, zMax = 3.0, zStep = 0.01, 
                 ZPError = 0.0, ZPOffsets = None, templatesDir = None, EBMinusVList = [0.0],
                 emLinesScaleList = [0.0]):
                 # EBMinusVList = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5],
                 # emLinesScaleList = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]):
                 # EBMinusVList = list(np.linspace(0, 0.48, 13))):
                 # EBMinusVList = list(np.linspace(0, 1.5, 16)),):
        """Sets up the stuff we would otherwise calculate every time, i.e., the templates.
        
        """

        if zStep < 0:
            print("... WARNING: zStep was negative - forced to be positive")
            zStep=abs(zStep)
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
            # DR10 adds i-band
            self.bands=['g', 'r', 'i', 'z', "w1", "w2"]
            for band in self.bands:
                if band[0] != "w":
                    self.passbandsList.append(astSED.Passband(passbandsDir+"decam.%s.am1p4.dat.txt" % (band)))
                else:
                    self.passbandsList.append(astSED.Passband(passbandsDir+"RSR-%s.EE.txt" % (band.upper()), inputUnits = "microns"))
        elif passbandSet == 'LSST':
            self.bands=['u', 'g', 'r', 'i', 'z', 'y']
            for band in self.bands:
                self.passbandsList.append(astSED.Passband(passbandsDir+band+"LSST.res"))
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
            if passbandSet == 'DECaLS':
                print(">>> Using default DECaLSDR10_learned_N20 templates [as DECaLS passbands used]")
                pickleFileName=None
                SEDDir=zCluster.__path__[0]+os.path.sep+'SED'
                self.SEDFiles=glob.glob(SEDDir+os.path.sep+"DECaLSDR10_learned_N20"+os.path.sep+"*.sed")
            else:
                # NOTE: This is the default template set used for every passbandSet in zCluster v0.5.1 and earlier
                print(">>> Using default EAZY_v1.0 and CWW templates")
                pickleFileName=None
                #SEDDir=resource_filename('zCluster', 'SED/')
                SEDDir=zCluster.__path__[0]+os.path.sep+'SED'
                self.SEDFiles=glob.glob(SEDDir+os.path.sep+"EAZY_v1.0"+os.path.sep+"*.dat")+ \
                            glob.glob(SEDDir+os.path.sep+"CWW"+os.path.sep+"*.sed")
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
            self.SEDFiles=self.SEDFiles+glob.glob(templatesDir+os.path.sep+"*spec.txt")
            if len(self.SEDFiles) == 0:
                raise Exception("No SED templates found - check path")
        # Because getting EOF errors under MPI
        pickleFileName=None
        if pickleFileName is not None and os.path.exists(pickleFileName):
            with open(pickleFileName, "rb") as pickleFile:
                unpickler=pickle.Unpickler(pickleFile)
                self.numModels=unpickler.load()
                self.modelFlux=unpickler.load()
                self.templateIndex=unpickler.load()
                self.modelFlux2=self.modelFlux**2
        else:
            self.numModels=len(self.SEDFiles)*len(EBMinusVList)*len(emLinesScaleList)
            i=0
            t=0
            self.templateIndex=[]
            for f in self.SEDFiles:
                s=astSED.SED()
                s.loadFromFile(f)
                t=t+1
                for emScaleFactor in emLinesScaleList:
                    if emScaleFactor > 0:
                        s.addEmissionLines(emScaleFactor)
                    for EBMinusV in EBMinusVList:
                        if EBMinusV > 0:
                            s.extinctionCalzetti(EBMinusV) # modifies z0flux, so this should be okay
                        for z in self.zRange:
                            s.redshift(z)
                            modelSEDDict=s.getSEDDict(self.passbandsList)
                            modelSEDDict['E(B-V)']=EBMinusV
                            modelSEDDict['ageGyr']=0.0
                            modelSEDDict['z']=z
                            modelSEDDict['fileName']=f
                            modelSEDDict['modelListIndex']=i
                            # modelSEDDict['SED']=s.copy() # Only add this if we need [e.g. rest frame colours] - uses lots of RAM
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
        if templatesDir is not None and pickleFileName is not None:
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
        #MStar=-21.5
        #MStarBand='r'
        #magsBrighterMStarCut=3  # Needs to be a little brighter to leave in BCGs
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
                
                try:
                    modFlux=np.average(modFlux, axis = 0, weights = chiSqProb)
                    modFluxJy=(modFlux*self.effLMicron2)/3e-13 
                    modMags=-2.5*np.log10(modFluxJy/1e23)-48.6
                    diffMags.append(magAB-modMags)
                except:
                    continue
            
        # Average after taking out differenced 99 values for missing data
        diffMags=np.array(diffMags)
        ZPOffsets=[]
        for i in range(diffMags.shape[1]):
            d=diffMags[:, i]
            ZPOffsets.append(np.median(d[np.less(d, 1)]))
        self.ZPOffsets=np.array(ZPOffsets)
        print("... offsets found: %s mag." % (str(self.ZPOffsets)))

    
    def calcPhotoRedshifts(self, galaxyCatalog, calcMLRedshiftAndOdds = False, returnPZ = True,
                           storeBestFitSEDModels = False):
        """Calculates photometric redshifts and adds to the galaxy catalog in place.
        
        NOTE: since normally we're normally only interested in p(z), this only returns
        the maximum likelihood z and BPZ-style odds parameter if calcMLRedshiftAndOdds = True.
        
        """
               
        # Note that we now aren't checking for junk/bad bands in the catalogue
        # That should have been dealt with before getting to this stage
        # Galaxies are fitted in batches using matrix products, rather than one at a time. For each galaxy,
        # with fluxes f, inverse variance weights w, and model fluxes M (one row per template/z combination):
        #   norm = sum(f*w*M)/sum(w*M^2)
        #   chiSq = sum(w*(f - norm*M)^2) = sum(w*f^2) - sum(f*w*M)^2/sum(w*M^2)
        # so both sums over bands can be done for a whole batch of galaxies in one matrix multiplication.
        if len(galaxyCatalog) == 0:
            return None
        ten23=(10**23.0)
        if self.ZPOffsets.sum() != 0:
            print("... applying zero-point offsets:", self.ZPOffsets, "...")

        # If we don't have a band, include a ridiculous mag with ridiculous error so zero weight in fit
        magAB=np.array([[galaxy[k] if k in galaxy else 99.0 for k in self.bands] for galaxy in galaxyCatalog],
                       dtype = float)
        magErrAB=np.array([[galaxy['%sErr' % (k)] if k in galaxy else 99.0 for k in self.bands]
                           for galaxy in galaxyCatalog], dtype = float)
        magAB=magAB-self.ZPOffsets

        # Incorporate a global photometric calibration error (applied to all bands)
        # Helps to get catalogs with tiny errors to behave
        magErrAB=np.sqrt(magErrAB**2+self.ZPError**2)

        fluxJy=ten23*np.power(10, (-(magAB+48.6)/2.5))
        sedFlux=3e-13*fluxJy/self.effLMicron2
        fluxJyErr=ten23*np.power(10, (-(magAB-magErrAB+48.6)/2.5))
        sedFluxErr=(3e-13*fluxJyErr/self.effLMicron2)-sedFlux
        with np.errstate(divide = 'ignore', invalid = 'ignore', over = 'ignore'):
            weights=1.0/sedFluxErr**2

        # Mag prior (one row per galaxy)
        absMag=magAB[:, self.magPriorBand].reshape(-1, 1)-5.0*np.log10(1e5*self.dlRange)
        if np.ndim(self.magPriorCut) == 0:
            pPrior=np.greater(absMag, self.magPriorCut)
        elif len(self.magPriorCut) == 2:
            pPrior=np.logical_and(np.greater(absMag, self.magPriorCut[0]), np.less(absMag, self.magPriorCut[1]))
        else:
            raise Exception("magPriorCut should have only 1 or 2 elements")

        # Batch size chosen to keep the (galaxies x models) work arrays to ~1 Mb each, so they stay in CPU cache
        # (in testing, this was ~3x faster than batches of ~40 Mb)
        numModelRows=self.modelFlux.shape[0]
        batchSize=max(1, int(1e5/numModelRows))
        modelFluxT=self.modelFlux.T
        modelFlux2T=self.modelFlux2.T
        numGals=len(galaxyCatalog)
        pzArray=np.zeros([numGals, self.zRange.shape[0]])
        for start in range(0, numGals, batchSize):
            end=min(start+batchSize, numGals)
            f=sedFlux[start:end]
            w=weights[start:end]
            with np.errstate(divide = 'ignore', invalid = 'ignore', over = 'ignore'):
                fw=f*w
                num=fw @ modelFluxT
                den=w @ modelFlux2T
                chiSq=np.sum(fw*f, axis = 1).reshape(-1, 1)-(num*num)/den
            chiSq[np.isnan(chiSq)]=1e6   # throw these out, should check this out and handle more gracefully

            # This is useful if we want to add rest frame color calculation
            if storeBestFitSEDModels == True:
                bestIndices=np.argmin(chiSq, axis = 1)
                for i, bestIndex in enumerate(bestIndices):
                    galaxyCatalog[start+i]['bestFitSED']=self.modelSEDDictList[self.templateIndex[bestIndex]]['SED']

            # This uses all templates at once - take max probability over templates at each z
            # NOTE: max(exp(-chiSq/2)) == exp(-min(chiSq)/2), so take the min first (far fewer exp calls)
            chiSq=chiSq.reshape([end-start, self.numModels, self.zRange.shape[0]])
            pzArray[start:end]=np.exp(-np.min(chiSq, axis = 1)/2)

        pzArray=pzArray*pPrior
        # Normalise
        pzNorm=np.trapezoid(pzArray, self.zRange, axis = 1)
        nonZero=pzNorm != 0
        pzArray[nonZero]=pzArray[nonZero]/pzNorm[nonZero].reshape(-1, 1)

        for galaxy, pz in zip(galaxyCatalog, pzArray):
            if calcMLRedshiftAndOdds == True:
                zp, odds=self.calculateMLRedshiftAndOdds(pz, dzOdds = 0.2)
                galaxy['zPhot']=zp
                galaxy['odds']=odds
            if returnPZ == True:
                galaxy['pz']=pz
                galaxy['pz_z']=self.zRange
        

    def calculateMLRedshiftAndOdds(self, pz, dzOdds = 0.2, method = 'max'):
        """Calculates maximum likelihood z and BPZ/EAZY style odds for given pz, zArray
        
        Method == 'max' finds z based on max peak in pz
        Method == 'odds' finds z based on z at which odds is maximised
        
        """
        
        # zRange is a uniform grid, so z -> index is linear (this replaces a spline fit made on every call)
        zGridStep=self.zRange[1]-self.zRange[0]
        def zToIndex(zs):
            return (np.array(zs)-self.zRange[0])/zGridStep
        z=self.zRange[np.argmax(pz)]

        if method == 'max':
            zMin=z-dzOdds
            zMax=z+dzOdds
            indexMin, indexMax=zToIndex([zMin, zMax])
            indexMin=int(round(indexMin))
            indexMax=int(round(indexMax))
            if indexMin < 0:
                indexMin=0
            if indexMax > pz.shape[0]-1:
                indexMax=pz.shape[0]-1
            odds=np.trapezoid(pz[indexMin:indexMax], self.zRange[indexMin:indexMax])
        
        elif method == 'odds':
            zOdds=[]
            for zi in self.zRange:
                zMin=zi-dzOdds
                zMax=zi+dzOdds
                indexMin, indexMax=zToIndex([zMin, zMax])
                indexMin=int(round(indexMin))
                indexMax=int(round(indexMax))
                if indexMin < 0:
                    indexMin=0
                if indexMax > pz.shape[0]-1:
                    indexMax=pz.shape[0]-1
                odds=np.trapezoid(pz[indexMin:indexMax], self.zRange[indexMin:indexMax])
                zOdds.append(odds)
            zOdds=np.array(zOdds)
            z=self.zRange[zOdds.tolist().index(zOdds.max())]
            odds=zOdds[zOdds.tolist().index(zOdds.max())]    
    
        return [z, odds]


    def estimateStellarMasses(self, galaxyCatalog, stellarMassModelDir, z = None):
        """Given a directory containing BC03-format stellar population models, estimate
        the stellar mass of galaxies in the given catalog.

        Args:
            galaxyCatalog (:obj:`list`): Galaxy catalog as a list of dictionaries, i.e.,
                in the format returned by self.calcPhotoRedshifts().
            stellarMassModelDir (:obj:`str`): Path to a directory containing the stellar
                population models (BC03 format for now).
            z (:obj:`float`, optional): If given, the redshift will be fixed to this value
                and applied to all the galaxies in the catalog (this is what you want for
                galaxy clusters, and is quicker). If None, then the maximum likelihood
                redshift of each individual galaxy will be used (this will be very slow,
                but actually isn't implemented yet...).

        Returns:
            None ['log10StellarMass' key is added in-place to each galaxy in galaxyCatalog]

        """

        if z is None:
            raise Exception("Stellar mass estimation at individual galaxy maximum likelihood redshifts is not implemented yet.")

        # Make model SEDs - this is very time consuming and needs a lot of speeding up...
        modelSEDDictList=sm.setUpStellarMassSEDs(stellarMassModelDir, self.passbandsList, z)

        # Fit each observed SED
        wantedKeys=['log10StellarMass']
        count=0
        DL=astCalc.dl(z)
        print(">>> Estimating stellar masses")
        for objDict in galaxyCatalog:
            count=count+1
            print("... %d/%d ..." % (count, len(galaxyCatalog)))
            mags=[]
            magErrs=[]
            for band in self.bands:
                if band in list(objDict.keys()):
                    mags.append(objDict[band])
                    magErrs.append(objDict[band+"Err"])
                else:
                    mags.append(99)
                    magErrs.append(99)
            obsSEDDict=astSED.mags2SEDDict(mags, magErrs, self.passbandsList)
            distNorm=4*np.pi*np.power(DL*3.08567758e24, 2)
            fitResult=sm.fitSEDDictAndCalcStellarMass(obsSEDDict, modelSEDDictList, distNorm)
            # Insert monte-carlo error estimation here...
            for key in wantedKeys:
                objDict[key]=fitResult[key]
