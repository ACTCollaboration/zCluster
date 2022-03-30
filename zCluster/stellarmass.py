"""

This module contains routines related to handling stellar population models,
for the purpose of stellar mass estimates.

"""

import os
import sys
import numpy as np
import astropy.table as atpy
import glob
import time
from astLib import *
from scipy import interpolate
from scipy import stats
import pylab
import itertools

#-------------------------------------------------------------------------------------------------------------
def setUpStellarMassSEDs(stellarMassModelDir, passbandsList, z, solarMetallicityOnly = True):
    """Set-up stellar population synthesis model SEDs for the purposes of stellar mass estimation.

    Only BC03 format is supported at the moment.

    """

    print(">>> Setting up stellar population models")
    t0=time.time()
    EBMinusVList=np.linspace(0, 0.48, 13).tolist()
    modelSEDDictList=[]
    if solarMetallicityOnly == True:
        modelSEDFileNames=glob.glob("%s/bc03_m62*.gplm" % (stellarMassModelDir))
    else:
        modelSEDFileNames=glob.glob("%s/bc03_*.gplm" % (stellarMassModelDir))
    count=0
    for modelSEDFileName in modelSEDFileNames:
        count=count+1
        print("... %d/%d ..." % (count, len(modelSEDFileNames)))
        modelSED=astSED.BC03Model(modelSEDFileName)
        # Only use ages less than age of the Universe
        ages=np.array(modelSED.ages)[np.less(modelSED.ages, astCalc.tz(z))]
        # Only load the mass file once per model SED - this is a bottle neck if we repeat for every age!
        modelMassDict=loadBC03MassFile(modelSEDFileName)
        #t00=time.time()
        for ageGyr in ages:
            s=modelSED.getSED(ageGyr, z = z)
            for EBMinusV in EBMinusVList:
                if EBMinusV > 0:
                    s.extinctionCalzetti(EBMinusV) # modifies z0flux, so this should be okay
                modelSEDDict=s.getSEDDict(passbandsList)
                modelSEDDict['labels']="label"
                modelSEDDict['EBMinusV']=EBMinusV
                modelSEDDict['ageGyr']=ageGyr
                modelSEDDict['z']=z
                modelSEDDict['fileName']=modelSED.fileName
                modelSEDDict['metallicity']=getMetallicity("BC03", modelSED.fileName)
                modelSEDDict['tauGyr']=getTauGyr("BC03", modelSED.fileName)
                modelSEDDict['modelListIndex']=count
                modelSEDDict['stellarMass']=interpolate.splev(modelSEDDict['ageGyr'], modelMassDict['tckMass'])
                modelSEDDictList.append(modelSEDDict)
    t1=time.time()
    print("... took %.3f sec" % (t1-t0))

    return modelSEDDictList

#-------------------------------------------------------------------------------------------------------------
def loadBC03MassFile(modelPath):
    """Loads BC03 *.4color files, which contain the stellar mass info, corresponding to the given model path.
    Returns dictionaries with numpy arrays of the different masses as a fn. of age., plus output from
    splrep which can be used to interpolate mass from age - e.g. interpolate.splev([ageGyr], dict['tckMass']

    Also extract tau from modelPath and use to calculate SFR evolution, returned in keys 'tauGyr', 'SFR'.
    Multiply SFR by stellar mass to get SFR in solar masses per year.

    Heck, while we're at it, parse the metallicity from modelPath into something sensible too.

    """

    massPath=modelPath.replace(".gplm", ".4color")
    inFile=open(massPath, "r")
    lines=inFile.readlines()
    inFile.close()
    age=[]
    Mstar=[]
    Mbol=[]
    MStarToLV=[]
    Vabs=[]    # abs mag in rest frame V through Buser V filter
    for line in lines:
        if line[0] != "#" and len(line) > 3:
            bits=line.split()
            if len(bits) > 3:   # incomplete lines for some CB07 files - not sure if this is a problem
                age.append((10.0**float(bits[0]))/1e9)
                Mstar.append(float(bits[6]))
                Mbol.append(float(bits[1]))
                MStarToLV.append(float(bits[5]))
                Vabs.append(float(bits[3]))
            #else:
                #print "What now?"
                #ipshell()
                #sys.exit()
    age=np.array(age)
    Mstar=np.array(Mstar)
    Mbol=np.array(Mbol)

    # Sanity check - does the mass file contain a similar range of ages to the model file?
    # If not, masses can be garbage as the interpolation gets extrapolated in a stupid fashion
    m=astSED.BC03Model(modelPath)
    if max(m.ages) > age.max():
        raise(Exception, "Mass file corresponding to %s is truncated (not enough ages)" % (modelPath))

    # All of the cubic interpolation here seems somewhat flakey for BC03/CB07?
    s=interpolate.splrep(age, Mstar, k=1)
    splMbol=interpolate.splrep(age, Mbol, k=1)
    splVabs=interpolate.splrep(age, Vabs, k=1)
    splMStarToLV=interpolate.splrep(age, MStarToLV, k=1)

    # Parse tau for tau models from file name according to my convention, set up interpolation from age too
    tauGyr=float(modelPath[modelPath.find("tau")+3:modelPath.find(".gplm")].replace("p", "."))
    SFR=(1.0/tauGyr)*np.exp(-age/tauGyr)/1e9
    sSFR=interpolate.splrep(age, SFR, k=1)

    # Parse metallicity into ZSun value - these from BC03 README file
    labels=["m22", "m32", "m42", "m52", "m62", "m72"]
    metallicities=[0.005, 0.02, 0.2, 0.4, 1.0, 2.5]
    labelInFileName=modelPath[modelPath.find("_m")+1:modelPath.find("_m")+4]
    ZSun=metallicities[labels.index(labelInFileName)]

    return {'ageGyr': age, 'Mstar': Mstar, 'ZSun': ZSun, 'tauGyr': tauGyr, 'tckMass': s, 'tckSFR': sSFR,
            'Mbol': Mbol, 'Vabs': Vabs, 'tckMbol': splMbol, 'tckVabs': splVabs, 'tckMStarToLV': splMStarToLV}

#-------------------------------------------------------------------------------------------------------------
def getTauGyr(modelFamily, modelPath):
    """Returns value of tau (SFH) by parsing the modelFileName

    """
    if modelFamily == "M05":
        tauGyr=float(modelPath[modelPath.find("_e")+2:modelPath.find("_z")])
    elif modelFamily == "BC03":
        tauGyr=float(modelPath[modelPath.find("tau")+3:modelPath.find(".gplm")].replace("p", "."))
    else:
        raise(Exception, "Need to add code to parse tau values from model file names for this family")

    return tauGyr

#-------------------------------------------------------------------------------------------------------------
def getMetallicity(modelFamily, modelPath):
    """Returns value of tau (SFH) by parsing the modelFileName

    """
    if modelFamily == "M05":
        print("no code to get M05 metallicity yet")
        ipshell()
        sys.exit()
        tauGyr=float(modelPath[modelPath.find("_e")+2:modelPath.find("_z")])
    elif modelFamily == "BC03":
        # From BC03 README
        ZSun=0.02
        ZMapDict={'m22': 0.0001/ZSun,
                  'm32': 0.0004/ZSun,
                  'm42': 0.004/ZSun,
                  'm52': 0.008/ZSun,
                  'm62': 0.02/ZSun,
                  'm72': 0.05/ZSun}
        fileNameKey=os.path.split(modelPath)[-1].split("_")[1]
        Z=ZMapDict[fileNameKey]
    else:
        raise(Exception, "Need to add code to parse tau values from model file names for this family")

    return Z

#------------------------------------------------------------------------------------------------------------
def fitSEDDictAndCalcStellarMass(SEDDict, modelSEDDictList, distNorm):
    """Custom version of the astSED routine that uses priors to throw out unphysical solutions

    Assuming BC03 here
    """

    modelFlux=[]
    modelMasses=[]
    for modelSEDDict in modelSEDDictList:
        modelFlux.append(modelSEDDict['flux'])
        modelMasses.append(modelSEDDict['stellarMass'])
    modelFlux=np.array(modelFlux)
    modelMasses=np.array(modelMasses)
    sedFlux=np.array([SEDDict['flux']]*len(modelFlux))
    sedFluxErr=np.array([SEDDict['fluxErr']]*len(modelFlux))

    # Analytic expression below is for normalisation at minimum chi squared (see note book)
    norm=np.sum((modelFlux*sedFlux)/(sedFluxErr**2), axis=1)/np.sum(modelFlux**2/sedFluxErr**2, axis=1)
    norms=np.array([norm]*modelFlux.shape[1]).transpose()
    chiSq=np.sum(((sedFlux-norms*modelFlux)**2)/sedFluxErr**2, axis=1)
    chiSq[np.isnan(chiSq)]=1e6   # throw these out, should check this out and handle more gracefully
    minChiSq=chiSq.min()
    bestMatchIndex=np.equal(chiSq, minChiSq).nonzero()[0][0]
    bestNorm=norm[bestMatchIndex]
    bestChiSq=minChiSq
    bestChiSqContrib=((sedFlux[bestMatchIndex]-norms[bestMatchIndex]*modelFlux[bestMatchIndex])**2)\
                        /sedFluxErr[bestMatchIndex]**2

    # Stellar mass
    stellarMasses=norm*distNorm*modelMasses
    bestStellarMass=stellarMasses[bestMatchIndex]

    fitResult={'minChiSq': bestChiSq,
                'chiSqContrib': bestChiSqContrib,
                'allChiSqValues': chiSq,
                'ageGyr': modelSEDDictList[bestMatchIndex]['ageGyr'],
                'norm': bestNorm,
                'log10StellarMass': np.log10(bestStellarMass),
                'z': modelSEDDictList[bestMatchIndex]['z'],
                'EBMinusV': modelSEDDictList[bestMatchIndex]['EBMinusV'],
                'tauGyr': modelSEDDictList[bestMatchIndex]['tauGyr'],
                'metallicity': modelSEDDictList[bestMatchIndex]['metallicity'],
                'bestFitModelSEDDict': modelSEDDictList[bestMatchIndex]}

    modelPath=modelSEDDictList[bestMatchIndex]['fileName']
    bestFitModel=astSED.BC03Model(modelPath)

    s=bestFitModel.getSED(fitResult['ageGyr'], z = fitResult['z'])
    s.extinctionCalzetti(fitResult['EBMinusV'])
    s.flux=s.flux*fitResult['norm']
    fitResult['bestFitSED']=s

    # Error bar by weighted chi-sq probability
    # NOTE: This seems to be the problem in crashing when multi-processing?
    #prob=stats.chisqprob(chiSq, 1)
    #logMasses=np.log10(stellarMasses)
    #massBinEdges=np.linspace(8.0, 13.5, 800)
    #massBinCentres=(massBinEdges[1]-massBinEdges[0])/2.0+massBinEdges[:-1]
    #binProb=np.zeros(len(massBinCentres))
    #for i in range(len(massBinCentres)):
        #binMin=massBinEdges[i]
        #binMax=massBinEdges[i+1]
        #mask=np.logical_and(np.greater(logMasses, binMin), np.less(logMasses, binMax))
        #binProb[i]=prob[mask].sum()
    #norm=np.trapz(binProb, massBinCentres)
    #binProb=binProb/norm
    #maxLikelihoodMass=massBinCentres[np.where(binProb == binProb.max())]
    #myProb=0.0
    #count=0
    #peakIndex=np.where(binProb == binProb.max())[0][0]
    #while myProb < 0.68:
        #count=count+1
        #myProb=np.trapz(binProb[peakIndex-count:peakIndex+count], massBinCentres[peakIndex-count:peakIndex+count])
    #minMass=massBinCentres[peakIndex-count]
    #maxMass=massBinCentres[peakIndex+count]
    #errLogStellarMass=np.max([maxMass-maxLikelihoodMass, maxLikelihoodMass-minMass])
    #fitResult['errLog10StellarMass_fromChiSq']=errLogStellarMass

    # If commented out above
    fitResult['errLog10StellarMass_fromChiSq']=-99

    return fitResult

#-------------------------------------------------------------------------------------------------------------
def makeSEDPlot(SEDDict, fitResult, plotFileName, title = ""):
    """Make a SED plot with the best fit shown

    """

    # Convert error bars to +/- format, so can avoid going -ve so matplotlib can plot them
    s=fitResult['bestFitSED']
    xMin=3500
    xMax=18000
    label="Age = %.1f Gyr, Z/Z$_{\odot}$ = %.3f, $\\tau$ = %.1f, $EBMinusV$ = %.2f, log($M^*$/M$_{\odot}$) = %.2f, $\chi^2$ = %.1f" % (fitResult['ageGyr'], fitResult['metallicity'], fitResult['tauGyr'], fitResult['EBMinusV'], fitResult['log10StellarMass'], fitResult['minChiSq'])
    pylab.plot(s.wavelength, s.flux, 'k-', label = label)
    errPlusMinus=np.zeros([2, SEDDict['fluxErr'].shape[0]])
    errPlusMinus[0]=SEDDict['fluxErr']
    errPlusMinus[1]=SEDDict['fluxErr']
    errMask=np.less(SEDDict['flux']-SEDDict['fluxErr'], 0)
    errPlusMinus[0][errMask]=SEDDict['flux'][errMask]*0.9999
    pylab.errorbar(SEDDict['wavelength'], SEDDict['flux'], yerr = errPlusMinus,
                    fmt = "ro", mec = "k", ecolor = "k", linewidth = 1.5, mew = 1.5, ms = 8)
    pylab.xlim(xMin, xMax)
    pylab.xlabel("Wavelength ($\AA{}$)")
    pylab.ylabel("Flux Density (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    pylab.legend(loc = 'upper left', prop = {'size': 9})
    rangeMask=np.logical_and(np.greater(s.wavelength, xMin), np.less(s.wavelength, xMax))
    yMin=s.flux[rangeMask].min()
    yMax=s.flux[rangeMask].max()*1.1
    pylab.ylim(yMin, yMax)
    pylab.title(title)
    pylab.savefig(plotFileName)
    pylab.close()

