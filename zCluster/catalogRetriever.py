"""
    Routines for retrieving photometric catalogs from SDSS, CFHTLenS etc. 
    and parsing them into a dictionary list format.

    By default, we cache the raw .csv output (or whatever) from each 
    database under $HOME/.zCluster/cache/, but the zCluster script can 
    also use a user-specifed cache location.

    Object dictionaries made by routines in here should have photometry 
    keys defined like 'u', 'uErr' etc., so that routines in 
    PhotoRedshiftEngine can understand them.
    
    ---

    Copyright 2016 Matt Hilton (matt.hilton@mykolab.com)
    
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
import astropy.table as atpy
import pyfits
import urllib
import urllib2
import pylab
import subprocess
import time
import zCluster
import requests 
from astropy.io.votable import parse_single_table 
import IPython

#-------------------------------------------------------------------------------------------------------------
CACHE_DIR=os.environ['HOME']+os.path.sep+".zCluster"+os.path.sep+"cache"

#-------------------------------------------------------------------------------------------------------------
def makeCacheDir():
    """Makes the cache dir where we store raw .csv output from the databases we query, if it doesn't already
    exist.
    
    """
    
    if os.path.exists(CACHE_DIR) == False:
        os.makedirs(CACHE_DIR)
    
#-------------------------------------------------------------------------------------------------------------
def S82Retriever(RADeg, decDeg, halfBoxSizeDeg = 20.2/60.0, optionsDict = {}):
    """Retrieves SDSS Stripe 82 photometry at the given position.
    
    """

    makeCacheDir()
    
    # We used PhotoObj for slit masks, because Galaxy misses some stuff, particularly at high-z
    #tableName="PhotoObj"
    tableName="Galaxy"
    
    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
                
    url = 'http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp'
            
    outFileName=cacheDir+os.path.sep+"S82_%.4f_%.4f_%.4f_%s.csv" % (RADeg, decDeg, halfBoxSizeDeg, tableName)
    print "... getting SDSS Stripe82 photometry (file: %s) ..." % (outFileName)
    
    # First, check if we previously downloaded a catalog and it went wrong - if it did, delete the file so
    # we can try fetching it again
    #if os.path.exists(outFileName) == True:
        #inFile=file(outFileName, "r")
        #lines=inFile.readlines()
        #inFile.close()
        #if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
            #os.remove(outFileName)
        
    if os.path.exists(outFileName) == False:
    
        print "... fetching from the internet ..."

        # Clean galaxy photometry query - note flags for r-band only, may want to change
        # Assuming for limits here on equator, so don't bother with cos(dec)
        #
        # We may want to add something that does multiple queries if we want a bigger area from the standard
        # SDSS query interface
        
        # Less conservative query - seems it's actually the galaxy view that's missing S82 objects
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        sql="""SELECT ra,dec,dered_u,dered_g,dered_r,dered_i,dered_z,Err_u,Err_g,Err_r,Err_i,Err_z,flags_r,run 
            FROM %s
            WHERE 
            (run=206 OR run=106) AND
            ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f 
            """ % (tableName, RAMin, RAMax, decMin, decMax)
        #print sql
        #AND Err_g < 0.3 AND Err_r < 0.3 AND Err_i < 0.3 

        # Old, conservative query that ends up missing a hell of a lot of galaxies that we can see in S82 i-band
        #sql="""SELECT ra,dec,dered_u,dered_g,dered_r,dered_i,dered_z,Err_u,Err_g,Err_r,Err_i,Err_z,flags_r,run 
            #FROM Galaxy 
            #WHERE 
            #(run=206 OR run=106) AND
            #ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f 
            #AND ((flags_r & 0x10000000) != 0) 
            #-- detected in BINNED1 
            #AND ((flags_r & 0x8100000c00a0) = 0) 
            #-- not NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP, SATURATED, 
            #-- or BAD_COUNTS_ERROR. 
            #-- if you want to accept objects with interpolation problems for PSF mags, 
            #-- change this to: AND ((flags_r & 0x800a0) = 0) 
            #AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2)) 
            #-- not DEBLEND_NOPEAK or small PSF error 
            #-- (substitute psfmagerr in other band as appropriate) 
            #AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0) 
            #-- not INTERP_CENTER or not COSMIC_RAY - omit this AND clause if you want to 
            #-- accept objects with interpolation problems for PSF mags.
            #AND Err_g < 0.2 AND Err_r < 0.2 AND Err_i < 0.2 
            #""" % (RAMin, RAMax, decMin, decMax)
            
        # Filter SQL so that it'll work
        fsql = ''
        for line in sql.split('\n'):
            fsql += line.split('--')[0] + ' ' + os.linesep;
        
        params=urllib.urlencode({'cmd': fsql, 'format': "csv"})
        response=None
        while response == None:
            try:
                response=urllib2.urlopen(url+'?%s' % (params))
            except:
                print "Network down? Waiting 30 sec..."
                time.sleep(30)
            
        lines=response.read()
        lines=lines.split("\n")

        outFile=file(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=file(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    
    # Parse .csv into catalog
    if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
        catalog=None
    elif lines[1] == 'ERROR\n':
        raise Exception, "Error returned by SDSS query - probably too many objects?"
    else:
        catalog=[]
        idCount=0
        for line in lines[1:]: # first line always heading
            if len(line) > 3:
                photDict={}
                idCount=idCount+1
                bits=line.replace("\n", "").split(",")
                photDict['id']=idCount    # just so we have something
                try:
                    photDict['RADeg']=float(bits[0])
                except:
                    if lines[1][:46] == '"ERROR: Maximum 60 queries allowed per minute.':
                        print "... exceeded server queries per minute limit - waiting ..."
                        time.sleep(70)
                        os.remove(outFileName)
                        return "retry"
                    else:
                        print "what?"
                        ipshell()
                        sys.exit()
                photDict['decDeg']=float(bits[1])
                photDict['u']=float(bits[2])
                photDict['g']=float(bits[3])
                photDict['r']=float(bits[4])
                photDict['i']=float(bits[5])
                photDict['z']=float(bits[6])
                photDict['uErr']=float(bits[7])
                photDict['gErr']=float(bits[8])
                photDict['rErr']=float(bits[9])
                photDict['iErr']=float(bits[10])
                photDict['zErr']=float(bits[11])
                photDict['run']=int(bits[13])
                # Apply mag error cuts if given
                # We're just making the mag unconstrained here (missing data), rather than applying a limit
                # If we don't have a minimum of three useful bands, reject
                if 'maxMagError' in optionsDict.keys():
                    keep=checkMagErrors(photDict, optionsDict['maxMagError'])
                else:
                    keep=True
                if keep == True:
                    catalog.append(photDict)
        
        # This guards against us letting through an empty catalog -  e.g. if not in Stripe 82, we get a 
        # different run number, so object doesn't get added to catalog.
        if catalog == []:
            catalog=None
            
    return catalog

#-------------------------------------------------------------------------------------------------------------
def checkMagErrors(photDict, maxMagError, minBands = 3, bands = ['u', 'g', 'r', 'i', 'z']):
    """Checks if the magnitudes in photDict are less than maxMagError, for a minimum number of minBands 
    photometric bands.
    
    Returns True if the photDict passes the test, False if not
    
    """
    
    keep=True
    rejectCount=0
    for band in bands:
        if photDict['%sErr' % (band)] > maxMagError or photDict['%s' % (band)] < 0:
            photDict['%s' % (band)]=99.
            photDict['%sErr' % (band)]=99.
            rejectCount=rejectCount+1
    if len(bands) - rejectCount < minBands:
        keep=False
    
    return keep

#-------------------------------------------------------------------------------------------------------------
def PS1Retriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves PS1 photometry at the given position.
    
    NOTE: not using y-band at the moment.
    
    """

    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)

    outFileName=cacheDir+os.path.sep+"PS1_%.4f_%.4f_%.4f.xml" % (RADeg, decDeg, halfBoxSizeDeg)      
    print "... getting PS1 photometry (file: %s) ..." % (outFileName)
    
    if os.path.exists(outFileName) == False or 'refetch' in optionsDict.keys() and optionsDict['refetch'] == True:
        print "... fetching from the internet ..."
        minDet=1
        r=requests.get('https://archive.stsci.edu/panstarrs/search.php', 
                        params= {'RA': RADeg, 'DEC': decDeg, 
                                 'SR': halfBoxSizeDeg, 'max_records': 30000, 
                                 'outputformat': 'VOTable',
                                 'ndetections': ('>%d' % minDet)}) 
        outFile=file(outFileName, 'w') 
        outFile.write(r.text) 
        outFile.close()
    
    tab=parse_single_table(outFileName).to_table(use_names_over_ids = True)

    # Parse table into catalog
    if len(tab) == 0:
        catalog=None
    else:
        EBMinusV=getEBMinusV(RADeg, decDeg, optionsDict = optionsDict) # assume same across field
        catalog=[]
        idCount=0
        # As per this page, use (rmeanpsfmag - rmeankronmag) >= 0.5 to remove stars
        # https://confluence.stsci.edu/display/PANSTARRS/PS1+Sample+queries#PS1Samplequeries-GalaxyCandidatesforK2C14SNSearch
        tab=tab[np.where(tab['rMeanPSFMag']-tab['rMeanKronMag'] >= 0.5)]
        for row in tab:
            idCount=idCount+1
            photDict={}
            photDict['id']=idCount    # just so we have something - we could PS1 ID but skipping for now
            photDict['RADeg']=row['raMean']
            photDict['decDeg']=row['decMean']
            photDict['g']=row['gMeanKronMag']
            photDict['r']=row['rMeanKronMag']
            photDict['i']=row['iMeanKronMag']
            photDict['z']=row['zMeanKronMag']
            photDict['gErr']=row['gMeanKronMagErr']
            photDict['rErr']=row['rMeanKronMagErr']
            photDict['iErr']=row['iMeanKronMagErr']
            photDict['zErr']=row['zMeanKronMagErr']
            
            # Correct for dust extinction
            photDict['g']=photDict['g']-EBMinusV*3.793 
            photDict['r']=photDict['r']-EBMinusV*2.751 
            photDict['i']=photDict['i']-EBMinusV*2.086 
            photDict['z']=photDict['z']-EBMinusV*1.479 
                
            # Apply mag error cuts if given
            # For PS1, missing values are -999 - our current checkMagErrors routine will fish those out
            # We're just making the mag unconstrained here (missing data), rather than applying a limit
            # If we don't have a minimum of three useful bands, reject
            if 'maxMagError' in optionsDict.keys():
                keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['g', 'r', 'i', 'z'])
            else:
                keep=True
            if keep == True:
                catalog.append(photDict)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, DR = 7, optionsDict = {}):
    """Retrieves SDSS main photometry at the given position.
    
    """
    
    # PhotoPrimary avoids star-galaxy classification problems at higher z, which we need when making slit masks
    # BUT contamination from stars etc. increases scatter in (z_spec - z_phot) by a fair bit (outliers)
    #tableName="PhotoPrimary"
    tableName="Galaxy"
    
    makeCacheDir()
    
    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    if DR == 7:
        url='http://cas.sdss.org/astrodr7/en/tools/search/x_sql.asp'
        outFileName=cacheDir+os.path.sep+"SDSSDR7_%.4f_%.4f_%.4f.csv" % (RADeg, decDeg, halfBoxSizeDeg)
        lineSkip=1
    elif DR == 8:
        url='http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp'
        outFileName=cacheDir+os.path.sep+"SDSSDR8_%.4f_%.4f_%.4f.csv" % (RADeg, decDeg, halfBoxSizeDeg)
        lineSkip=1
    elif DR == 10:      
        url='http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'
        outFileName=cacheDir+os.path.sep+"SDSSDR10_%.4f_%.4f_%.4f.csv" % (RADeg, decDeg, halfBoxSizeDeg)
        lineSkip=2
    elif DR == 12:
        # For some reason, SDSS changed their whole web API in ~May 2016 without calling it a new DR
        #url='http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx'        
        #url='http://skyserver.sdss.org/dr12/en/tools/search/x_results.aspx'
        url='http://skyserver.sdss.org/dr12/en/tools/search/x_results.aspx?searchtool=SQL&TaskName=Skyserver.Search.SQL&syntax=NoSyntax&ReturnHtml=false&'
        outFileName=cacheDir+os.path.sep+"SDSSDR12_%.4f_%.4f_%.4f.csv" % (RADeg, decDeg, halfBoxSizeDeg)
        lineSkip=2		
   
    outFileName=outFileName.replace(".csv", "_%s.csv" % (tableName))
                                    
    print "... getting SDSS DR%d photometry (file: %s) ..." % (DR, outFileName)

    # First, check if we previously downloaded a catalog and it went wrong - if it did, delete the file so
    # we can try fetching it again
    #if os.path.exists(outFileName) == True:
        #inFile=file(outFileName, "r")
        #lines=inFile.readlines()
        #inFile.close()
        #if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
            #os.remove(outFileName)
        
    if os.path.exists(outFileName) == False or 'refetch' in optionsDict.keys() and optionsDict['refetch'] == True:
        
        print "... fetching from the internet ..."
        
        # Clean galaxy photometry query - note flags for r-band only, may want to change
        # We may want to add something that does multiple queries if we want a bigger area from the standard
        # SDSS query interface
        # NOTE: disabled the 'clean' photometry flags on 21/09/2015
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        # Sanity check: if RAMin, RAMax have order reversed, SDSS gives us nothing
        if RAMin > RAMax:
            print "RAMin, RAMax reversed"
            IPython.embed()
            sys.exit()
            newRAMax=RAMin
            RAMin=RAMax
            RAMax=newRAMax
        if decMin > decMax:
            print "decMin, decMax reversed"
            IPython.embed()
            sys.exit()
            newDecMax=decMin
            decMin=decMax
            decMax=newDecMax
        # Was FROM Galaxy,and taken out mag err constrain on 23 Oct 2015
        # PhotoPrimary is like PhotoObj but without multiple matches
        sql="""SELECT ra,dec,dered_u,dered_g,dered_r,dered_i,dered_z,Err_u,Err_g,Err_r,Err_i,Err_z,flags_r,run 
            FROM %s
            WHERE 
            ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f 
            -- AND ((flags_r & 0x10000000) != 0) 
            -- detected in BINNED1 
            -- AND ((flags_r & 0x8100000c00a0) = 0) 
            -- not NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP, SATURATED, 
            -- or BAD_COUNTS_ERROR. 
            -- if you want to accept objects with interpolation problems for PSF mags, 
            -- change this to: AND ((flags_r & 0x800a0) = 0) 
            -- AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2)) 
            -- not DEBLEND_NOPEAK or small PSF error 
            -- (substitute psfmagerr in other band as appropriate) 
            -- AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0) 
            -- not INTERP_CENTER or not COSMIC_RAY - omit this AND clause if you want to 
            -- accept objects with interpolation problems for PSF mags.
            --AND Err_g < 0.5 AND Err_r < 0.5 AND Err_i < 0.5 
            """ % (tableName, RAMin, RAMax, decMin, decMax)
    
        # Filter SQL so that it'll work
        fsql = ''
        for line in sql.split('\n'):
            fsql += line.split('--')[0] + ' ' + os.linesep;
    
        params=urllib.urlencode({'cmd': fsql, 'format': "csv"})
        response=None
        while response == None:
            try:
                response=urllib2.urlopen(url+'%s' % (params))
            except:
                print "Network down? Waiting 30 sec..."
                time.sleep(30)

        lines=response.read()
        lines=lines.split("\n")
        
        outFile=file(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=file(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    
    # Parse .csv into catalog
    if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
        catalog=None
    else:
        catalog=[]
        idCount=0
        for line in lines[lineSkip:]: # first line or two (depending on DR) always heading
            if len(line) > 3:
                photDict={}
                idCount=idCount+1
                bits=line.replace("\n", "").split(",")
                photDict['id']=idCount    # just so we have something
                try:
                    photDict['RADeg']=float(bits[0])
                except:
                    if lines[1][:46] == '"ERROR: Maximum 60 queries allowed per minute.':
                        print "... exceeded server queries per minute limit - waiting ..."
                        time.sleep(70)
                        os.remove(outFileName)
                        return "retry"
                    else:
                        if line.find("<html>") != -1:
                            time.sleep(70)
                            os.remove(outFileName)
                            return "retry"
                        else:
                            print "what?"
                            IPython.embed()
                            sys.exit()
                photDict['decDeg']=float(bits[1])
                photDict['u']=float(bits[2])
                photDict['g']=float(bits[3])
                photDict['r']=float(bits[4])
                photDict['i']=float(bits[5])
                photDict['z']=float(bits[6])
                photDict['uErr']=float(bits[7])
                photDict['gErr']=float(bits[8])
                photDict['rErr']=float(bits[9])
                photDict['iErr']=float(bits[10])
                photDict['zErr']=float(bits[11])
                photDict['run']=int(bits[13])
                # Apply mag error cuts if given
                # We're just making the mag unconstrained here (missing data), rather than applying a limit
                # If we don't have a minimum of three useful bands, reject
                if 'maxMagError' in optionsDict.keys():
                    keep=checkMagErrors(photDict, optionsDict['maxMagError'])
                else:
                    keep=True
                if keep == True:
                    catalog.append(photDict)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def DECaLSRetriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves DECaLS DR3 tractor catalogs (if they exist) at the given position
    
    """

    makeCacheDir()
    
    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)

    # Organised such that after this, have subdir with degrees RA (e.g. 000/ is 0 < RADeg < 1 degree)
    basePath="http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr3/tractor/"

    outFileName=cacheDir+os.path.sep+"AllDR3TractorCats_%.6f_%.6f_%.2f.fits" % (RADeg, decDeg,
                                                                            halfBoxSizeDeg)

    print "... getting DECaLS photometry (file: %s) ..." % (outFileName)

    if os.path.exists(outFileName) == False:
        
        bricksTab=atpy.Table(zCluster.__path__[0]+os.path.sep+"data"+os.path.sep+"decals-bricks.fits")

        # Find matching tractor catalogs and download (temporary) .fits table files
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        mask=np.logical_and(np.greater(bricksTab['ra1'], RAMin), np.less(bricksTab['ra2'], RAMax))
        mask=np.logical_and(mask, np.greater(bricksTab['dec1'], decMin))
        mask=np.logical_and(mask, np.less(bricksTab['dec2'], decMax))
        matchTab=bricksTab.where(mask)
        count=0
        tractorTabs=[]
        cacheFileNames=[]
        for mrow in matchTab:
            count=count+1
            if mrow['has_catalog'] == 1:
                print "... retrieving tractor catalog from web ..."
                url=basePath+"%03d" % np.floor(RADeg)
                url=url+os.path.sep+"tractor-"+mrow['brickname']+".fits"
                fileName=cacheDir+os.path.sep+"tractor-tmp_%d.fits" % (count)
                cacheFileNames.append(fileName)
                urllib.urlretrieve(url, filename = fileName)
                try:
                    tractorTab=atpy.Table(fileName)
                    tractorTabs.append(tractorTab)
                except:
                    print "... probably a 404 error for %s ..." % (fileName)
        # Stitch catalogs together and write to cache, delete temporary files
        if len(tractorTabs) > 0:
            uberTab=atpy.Table()
            for key in tractorTabs[0].keys():
                # Skip DECam aperture fluxes for now (8 x 6 arrays)
                #if key.find("decam_apflux") == -1:
                arr=[]
                for tracTab in tractorTabs:
                    for row in tracTab:
                        # skip over, e.g., decam_apflux
                        if type(row[key]) == np.ndarray and row[key].ndim == 2:
                            break
                        arr.append(row[key])
                if len(arr) > 0:
                    uberTab.add_column(key, arr)
            uberTab.table_name="DECaLS"
            uberTab.write(outFileName)
            for fileName in cacheFileNames:
                os.remove(fileName)
    
    # Load catalog
    if os.path.exists(outFileName) == True:
        
        tab=atpy.Table(outFileName)
 
        # Remove stars
        tab=tab.where(tab['type'] != 'PSF')

        # We could use DECaLS extinction coeffs, but for now do the same as we are elsewhere   
        # Skipping Y for now
        corrDict={'u': 5.155, 'g': 3.793, 'r': 2.751, 'i': 2.086, 'z': 1.479, 'Ks': 0.367}    
        EBMinusV=getEBMinusV(RADeg, decDeg, optionsDict = optionsDict) # assume same across field (much faster)

        catalog=[]
        for row in tab: 
            photDict={}
            photDict['id']=row['objid']
            photDict['RADeg']=row['ra']
            photDict['decDeg']=row['dec']
            # Flux in DECaLS is in nanomaggies: an object with AB mag = 22.5 has flux 1 nanomaggie
            # So, mag AB = -2.5 log (nanomaggies) + 22.5
            # Skipping Y for now
            bands=['u', 'g', 'r', 'i', 'z']# , 'Y']
            for b in bands:
                dustCorrMag=EBMinusV*corrDict[b]
                nanoMaggies=row['decam_flux'][bands.index(b)]
                if nanoMaggies <= 0:
                    mag=99.
                    magErr=99.
                else:
                    mag=-2.5*np.log10(nanoMaggies)+22.5
                    fluxErr=np.sqrt(1./row['decam_flux_ivar'][bands.index(b)])
                    magErr=1./(nanoMaggies/fluxErr)
                photDict[b]=mag-dustCorrMag
                photDict[b+"Err"]=magErr

            if 'maxMagError' in optionsDict.keys():
                keep=checkMagErrors(photDict, optionsDict['maxMagError'])
            else:
                keep=True
            if keep == True:
                catalog.append(photDict)
    
    else:
        catalog=None
                
    return catalog

#-------------------------------------------------------------------------------------------------------------
def SDSSDR7Retriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, optionsDict = {}):
    """Retrieves SDSS DR7 main photometry at the given position.
    
    """
    
    stuff=SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = halfBoxSizeDeg, DR = 7, optionsDict = optionsDict)
    return stuff

#-------------------------------------------------------------------------------------------------------------
def SDSSDR8Retriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, optionsDict = {}):
    """Retrieves SDSS DR8 photometry at the given position.
    
    """
    
    stuff=SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = halfBoxSizeDeg, DR = 8, optionsDict = optionsDict)
    return stuff

#-------------------------------------------------------------------------------------------------------------
def SDSSDR10Retriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, optionsDict = {}):
    """Retrieves SDSS DR10 photometry at the given position.
    
    """
    
    stuff=SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = halfBoxSizeDeg, DR = 10, optionsDict = optionsDict)
    return stuff

#-------------------------------------------------------------------------------------------------------------
def SDSSDR12Retriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves SDSS DR12 photometry at the given position.
    
    """
    
    stuff=SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = halfBoxSizeDeg, DR = 12, optionsDict = optionsDict)
    return stuff

#-------------------------------------------------------------------------------------------------------------
def CFHTLenSRetriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves CFHTLenS photometry, which works differently to other CFHT catalogues. 
    
    halfBoxSizeDeg is actually a radius in this case.
    
    NOTE: optionsDict must include 'maxMagError' key.
    
    """
    
    makeCacheDir()
    
    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    # CFHTLens - note: mag 99s where missing, this does star/galaxy sep based on lensfit fitclass
    #url="http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLens/cgi/queryt.pl?REQUEST=doQuery&LANG=ADQL&method=sync&format=ascii&query=SELECT%0D%0Aid%2C+ALPHA_J2000%2C+DELTA_J2000%2C+fitclass%2C+MAG_u%2C+MAGERR_u%2C+MAG_g%2C+MAGERR_g%2C+MAG_r%2C+MAGERR_r%2C+MAG_i%2C+MAGERR_i%2C+MAG_y%2C+MAGERR_y%2C+MAG_z%2C+MAGERR_z%0D%0AFROM%0D%0Acfht.clens%0D%0AWHERE%0D%0Afitclass%3E%3D0%0D%0AAND+fitclass%3C%3D0%0D%0AAND+contains%28pos%2Ccircle%28%27ICRS+GEOCENTER%27%2C$RADEG%2C$DECDEG%2C$SEARCHRAD%29%29%3D1%0D%0A" 
    # No star-galaxy separation
    #http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/community/CFHTLens/queryt.pl?REQUEST=doQuery&LANG=ADQL&method=sync&format=ascii&query=SELECT%0D%0Atop+10%0D%0Aid%0D%0AFROM%0D%0Acfht.clens%0D%0A
    # Updated March 2016: changed URLs at CADC
    url="http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/community/CFHTLens/queryt.pl?REQUEST=doQuery&LANG=ADQL&method=sync&format=ascii&query=SELECT%0D%0Aid%2C+ALPHA_J2000%2C+DELTA_J2000%2C+fitclass%2C+MAG_u%2C+MAGERR_u%2C+MAG_g%2C+MAGERR_g%2C+MAG_r%2C+MAGERR_r%2C+MAG_i%2C+MAGERR_i%2C+MAG_y%2C+MAGERR_y%2C+MAG_z%2C+MAGERR_z%0D%0AFROM%0D%0Acfht.clens%0D%0AWHERE%0D%0A+contains%28pos%2Ccircle%28%27ICRS+GEOCENTER%27%2C$RADEG%2C$DECDEG%2C$SEARCHRAD%29%29%3D1%0D%0A" 
    url=url.replace("$RADEG", "%.6f" % (RADeg))
    url=url.replace("$DECDEG", "%.6f" % (decDeg))
    url=url.replace("$SEARCHRAD", "%.6f" % (halfBoxSizeDeg))

    outFileName=cacheDir+os.path.sep+"CFHTLenS_%.4f_%.4f_%.4f.csv" % (RADeg, decDeg, halfBoxSizeDeg)

    print "... getting CFHTLenS photometry (file: %s) ..." % (outFileName) 

    # First, check if we previously downloaded a catalog and it went wrong - if it did, delete the file so
    # we can try fetching it again
    #if os.path.exists(outFileName) == True:
        #inFile=file(outFileName, "r")
        #lines=inFile.readlines()
        #inFile.close()
        #if len(lines) < 10:
            #os.remove(outFileName)
        
    if os.path.exists(outFileName) == False:
        
        print "... fetching from the internet ..."
            
        response=None
        while response == None:
            try:
                response=urllib2.urlopen(url)
            except:
                print "Network down, or CADC changed URLs again? Waiting 30 sec..."
                time.sleep(30)
                
        lines=response.read()
        lines=lines.split("\n")

        outFile=file(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=file(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()

    # Parse .csv into catalog
    if len(lines) < 10:
        catalog=None
    else:
        catalog=[]
        idCount=0
        #EBMinusV=getEBMinusV(RADeg, decDeg) # Extinction already applied in CFHTLenS mag columns
        for line in lines[1:]: # first line always heading
            if len(line) > 3:
                photDict={}
                idCount=idCount+1
                bits=line.replace("\n", "").split()                
                photDict['id']=idCount    # just so we have something
                try:
                    photDict['RADeg']=float(bits[1])
                except:
                    print "what?"
                    IPython.embed()
                    sys.exit()
                # Note: skip mag y, old version of i-band filter
                photDict['decDeg']=float(bits[2])
                photDict['u']=float(bits[4])
                photDict['g']=float(bits[6])
                photDict['r']=float(bits[8])
                photDict['i']=float(bits[10])
                photDict['z']=float(bits[14])
                photDict['uErr']=float(bits[5])
                photDict['gErr']=float(bits[7])
                photDict['rErr']=float(bits[9])
                photDict['iErr']=float(bits[11])
                photDict['zErr']=float(bits[15])
                                
                # Correct for dust extinction - not necessary for CFHTLenS, already applied
                #photDict['u']=photDict['u']-EBMinusV*5.155 
                #photDict['g']=photDict['g']-EBMinusV*3.793 
                #photDict['r']=photDict['r']-EBMinusV*2.751 
                #photDict['i']=photDict['i']-EBMinusV*2.086 
                #photDict['z']=photDict['z']-EBMinusV*1.479 
                
                # Apply mag error cuts if given
                # We're just making the mag unconstrained here (missing data), rather than applying a limit
                # If we don't have a minimum of three useful bands, reject
                if 'maxMagError' in optionsDict.keys():
                    keep=checkMagErrors(photDict, optionsDict['maxMagError'])
                else:
                    keep=True
                # If we want to keep what CFHTLenS classify as galaxies (galaxies are fitClass == 0)
                # NOTE: this makes no significant difference in our testing, so disabled
                #fitClass=int(bits[3])
                #if fitClass == 0:
                    #keep=True
                #else:
                    #keep=False
                if keep == True:
                    catalog.append(photDict)
                    
    return catalog

#-------------------------------------------------------------------------------------------------------------
def CFHTRetriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, survey = 'deep', optionsDict = {}):
    """Retrieves CFHT photometry  fields (all public now) using Steven Gwyn's CADC pages. 
    
    halfBoxSizeDeg is actually a radius in this case.
    
    """
     
    makeCacheDir()
    
    # Below describes the catalogs - they are done in matched apertures but are i-band selected
    #http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/docs/cfhtlsproc.html#photcat
    
    #url='http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/catquery.pl?method=direct&survey=deep&filter=i&format=ascii&select=+RA,+DEC,+U_MAG_AUTO,+U_MAGERR_AUTO,+G_MAG_AUTO,+G_MAGERR_AUTO,+R_MAG_AUTO,+R_MAGERR_AUTO,+I_MAG_AUTO,+I_MAGERR_AUTO,+Z_MAG_AUTO,+Z_MAGERR_AUTO&where=incirc($RA,$DEC,$SEARCHRAD)%0D%0AAND+i_MAG_AUTO+%3C%3D+i_maglimit%0D%0AAND+i_dubious+%3D0&bright=no&masked=no&ofield=all&col=U_MAG_AUTO&col=U_MAGERR_AUTO&col=G_MAG_AUTO&col=G_MAGERR_AUTO&col=R_MAG_AUTO&col=R_MAGERR_AUTO&col=I_MAG_AUTO&col=I_MAGERR_AUTO&col=Z_MAG_AUTO&col=Z_MAGERR_AUTO'
    # below is attempt at cleaner photometry
    #url='http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/catquery.pl?method=direct&survey=$SURVEY&filter=i&format=ascii&select=+RA,+DEC,+U_MAG_AUTO,+U_MAGERR_AUTO,+G_MAG_AUTO,+G_MAGERR_AUTO,+R_MAG_AUTO,+R_MAGERR_AUTO,+I_MAG_AUTO,+I_MAGERR_AUTO,+Z_MAG_AUTO,+Z_MAGERR_AUTO&where=incirc($RA,$DEC,$SEARCHRAD)%0D%0AAND+i_MAG_AUTO+%3C%3D+i_maglimit%0D%0AAND+i_dubious+%3D0%0D%0AAND+r_dubious+%3D0&%0D%0AAND+g_dubious+%3D0%0D%0AAND+u_dubious+%3D0%0D%0AAND+z_dubious+%3D0bright=no&masked=no&ofield=all&col=U_MAG_AUTO&col=U_MAGERR_AUTO&col=G_MAG_AUTO&col=G_MAGERR_AUTO&col=R_MAG_AUTO&col=R_MAGERR_AUTO&col=I_MAG_AUTO&col=I_MAGERR_AUTO&col=Z_MAG_AUTO&col=Z_MAGERR_AUTO'
    
    # March 2016: new query location
    #http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/megapipe/queryt.pl?slang=en&REQUEST=doQuery&LANG=ADQL&meth=sync&format=ascii&query=SELECT%0D%0ATOP+10%0D%0A+++cfhtlsID%0D%0AFROM%0D%0A+++cfht.deepi%0D%0AWHERE%0D%0A+++i_MAG_AUTO+%3C%3D+i_maglimit%0D%0A+++AND+i_dubious+%3D0%0D%0A

    print "update query"
    IPython.embed()
    sys.exit()
    
    # below should be cleaner still (5 sigma in g, r, i), plus actually have star-galaxy separation in
    url='http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/catquery.pl?method=direct&survey=$SURVEY&filter=i&format=ascii&select=+RA,+DEC,+U_MAG_AUTO,+U_MAGERR_AUTO,+G_MAG_AUTO,+G_MAGERR_AUTO,+R_MAG_AUTO,+R_MAGERR_AUTO,+I_MAG_AUTO,+I_MAGERR_AUTO,+Z_MAG_AUTO,+Z_MAGERR_AUTO,+I_CLASS_STAR&where=incirc($RA,$DEC,$SEARCHRAD)%0D%0AAND+i_MAG_AUTO+%3C%3D+i_maglimit%0D%0AAND+i_MAGERR_AUTO+%3C+0.2%0D%0AAND+r_MAGERR_AUTO+%3C+0.2%0D%0AAND+g_MAGERR_AUTO+%3C+0.2%0D%0AAND+I_CLASS_STAR+%3C+059%0D%0AAND+i_dubious+%3D0%0D%0AAND+r_dubious+%3D0&%0D%0AAND+g_dubious+%3D0%0D%0AAND+u_dubious+%3D0%0D%0AAND+z_dubious+%3D0bright=no&masked=no&ofield=all&col=U_MAG_AUTO&col=U_MAGERR_AUTO&col=G_MAG_AUTO&col=G_MAGERR_AUTO&col=R_MAG_AUTO&col=R_MAGERR_AUTO&col=I_MAG_AUTO&col=I_MAGERR_AUTO&col=Z_MAG_AUTO&col=Z_MAGERR_AUTO'
    
    # Aper mags don't help at all
    #url='http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/catquery.pl?method=direct&survey=$SURVEY&filter=i&format=ascii&select=+RA,+DEC,+U_MAG_APER_2,+U_MAGERR_APER_2,+G_MAG_APER_2,+G_MAGERR_APER_2,+R_MAG_APER_2,+R_MAGERR_APER_2,+I_MAG_APER_2,+I_MAGERR_APER_2,+Z_MAG_APER_2,+Z_MAGERR_APER_2&where=incirc($RA,$DEC,$SEARCHRAD)%0D%0AAND+i_MAG_AUTO+%3C%3D+i_maglimit%0D%0AAND+i_MAGERR_AUTO+%3C+0.2%0D%0AAND+r_MAGERR_AUTO+%3C+0.2%0D%0AAND+g_MAGERR_AUTO+%3C+0.2%0D%0AAND+R_CLASS_STAR+%3C+0.9%0D%0AAND+i_dubious+%3D0%0D%0AAND+r_dubious+%3D0&%0D%0AAND+g_dubious+%3D0%0D%0AAND+u_dubious+%3D0%0D%0AAND+z_dubious+%3D0bright=no&masked=no&ofield=all&col=U_MAG_AUTO&col=U_MAGERR_AUTO&col=G_MAG_AUTO&col=G_MAGERR_AUTO&col=R_MAG_AUTO&col=R_MAGERR_AUTO&col=I_MAG_AUTO&col=I_MAGERR_AUTO&col=Z_MAG_AUTO&col=Z_MAGERR_AUTO'
    
    # CFHTLens - note: mag 99s where missing, this does star/galaxy sep based on lensfit fitclass
    #http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLens/cgi/queryt.pl?REQUEST=doQuery&LANG=ADQL&method=sync&format=ascii&query=SELECT%0D%0Aid%2C+ALPHA_J2000%2C+DELTA_J2000%2C+MAG_u%2C+MAGERR_u%2C+MAG_g%2C+MAGERR_g%2C+MAG_r%2C+MAGERR_r%2C+MAG_i%2C+MAGERR_i%2C+MAG_y%2C+MAGERR_y%2C+MAG_z%2C+MAGERR_z%0D%0AFROM%0D%0Acfht.clens%0D%0AWHERE%0D%0Acontains%28pos%2Ccircle%28%27ICRS+GEOCENTER%27%2C36.1333%2C-7.5048%2C0.1500%29%29%3D1+and+fitclass%3D0%0D%0A
    
    url=url.replace("$SURVEY", survey)
    url=url.replace("$RA", "%.6f" % (RADeg))
    url=url.replace("$DEC", "%.6f" % (decDeg))
    url=url.replace("$SEARCHRAD", "%.6f" % (halfBoxSizeDeg))
    
    if survey == 'deep':
        label='CFHTDeep'
    elif survey == 'wide':
        label='CFHTWide'
    outFileName=CACHE_DIR+os.path.sep+"%s_%.4f_%.4f_%.4f.csv" % (label, RADeg, decDeg, halfBoxSizeDeg)

    print "... getting CFHTLS %s photometry (file: %s) ..." % (survey, outFileName) 

    # First, check if we previously downloaded a catalog and it went wrong - if it did, delete the file so
    # we can try fetching it again
    #if os.path.exists(outFileName) == True:
        #inFile=file(outFileName, "r")
        #lines=inFile.readlines()
        #inFile.close()
        #if len(lines) < 10:
            #os.remove(outFileName)
        
    if os.path.exists(outFileName) == False:
            
        response=urllib2.urlopen(url)
        lines=response.read()
        lines=lines.split("\n")

        outFile=file(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=file(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()

    # Parse .csv into catalog
    if len(lines) < 10:
        catalog=None
    else:
        catalog=[]
        idCount=0
        EBMinusV=getEBMinusV(RADeg, decDeg, optionsDict = optionsDict) # assume same across field (much faster)
        for line in lines[1:]: # first line always heading
            if len(line) > 3:
                photDict={}
                idCount=idCount+1
                bits=line.replace("\n", "").split()
                photDict['id']=idCount    # just so we have something
                try:
                    photDict['RADeg']=float(bits[0])
                except:
                    print "what?"
                    IPython.embed()
                    sys.exit()
                photDict['decDeg']=float(bits[1])
                photDict['u']=float(bits[2])
                photDict['g']=float(bits[4])
                photDict['r']=float(bits[6])
                photDict['i']=float(bits[8])
                photDict['z']=float(bits[10])
                photDict['uErr']=float(bits[3])
                photDict['gErr']=float(bits[5])
                photDict['rErr']=float(bits[7])
                photDict['iErr']=float(bits[9])
                photDict['zErr']=float(bits[11])
                photDict['classStar']=float(bits[12]) # 
                
                # Correct for dust extinction
                photDict['u']=photDict['u']-EBMinusV*5.155 
                photDict['g']=photDict['g']-EBMinusV*3.793 
                photDict['r']=photDict['r']-EBMinusV*2.751 
                photDict['i']=photDict['i']-EBMinusV*2.086 
                photDict['z']=photDict['z']-EBMinusV*1.479 
                
                catalog.append(photDict)
                
    # If we want to test seeing what our catalogs look like
    #ri=[]
    #gr=[]
    #for obj in catalog:
        #ri.append(obj['r']-obj['i'])
        #gr.append(obj['g']-obj['r'])
    #pylab.plot(gr, ri, '.')
    #pylab.xlim(-0.5, 1.5)
    #pylab.ylim(-0.5, 3)
    #ipshell()
    #sys.exit()
    
    # If we also wanted to retrieve image
    # Retrieve 5' i-band image of 2215 from CFHTLS - this gives a page with link to image to download
    #http://www4.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/cutcfhtls.pl?ra=333.99758911&dec=-17.6376571655%20&object=D3&size=5&units=arcminutes&F=I&wide=true&deep=true
    #ipshell()
    #sys.exit()
    
    return catalog
    
#-------------------------------------------------------------------------------------------------------------
def CFHTDeepRetriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, optionsDict = {}):
    """Retrieves CFHT deep photometry  fields (all public now) using Steven Gwyn's CADC pages. 
    
    halfBoxSizeDeg is actually a radius in this case.
    
    """ 
    
    catalog=CFHTRetriever(RADeg, decDeg, halfBoxSizeDeg, survey = 'deep')
    return catalog

#-------------------------------------------------------------------------------------------------------------
def CFHTWideRetriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, optionsDict = {}):
    """Retrieves CFHT wide photometry  fields (all public now) using Steven Gwyn's CADC pages. 
    
    halfBoxSizeDeg is actually a radius in this case.
    
    """
    
    catalog=CFHTRetriever(RADeg, decDeg, halfBoxSizeDeg, survey = 'wide')
    return catalog

#-------------------------------------------------------------------------------------------------------------
def getEBMinusV(RADeg, decDeg, optionsDict = {}):
    """Get E(B-V) due to Galactic dust from Schlegel maps, using IRSA web service
    
    """

    makeCacheDir()
    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
        
    fileName=cacheDir+os.path.sep+"SchlegelIRSA_%.6f_%.6f.xml" % (RADeg, decDeg)
    if os.path.exists(fileName) == False:
        url="http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr="+str(RADeg)+"%20"+str(decDeg)+"%20Equ%20J2000"
        urllib.urlretrieve(url, filename = fileName)
    
    inFile=file(fileName, "r")
    lines=inFile.readlines()
    inFile.close()
    for i in range(len(lines)):
        if lines[i].find("refPixelValueSFD") != -1:
            break
    EBMinusV=float(lines[i+1].split(" (")[0])

    return EBMinusV

#-------------------------------------------------------------------------------------------------------------
def parseFITSPhotoTable(tab, fieldIDKey = None, optionsDict = {}):
    """Parse .fits table into catalog list of dictionaries format. Extinction correction due to Galactic
    dust will be applied, using Schlegel dust maps. If fieldIDKey == None (the default), this is done
    at the mean RA, dec coords. Otherwise, sources are grouped by fieldID and the mean correction for
    each field is applied.
     
    Returns catalog (list of dictionaries)
    
    """

    # Options we may wish to play with
    magKey="MAG_AUTO"
    magErrKey="MAGERR_AUTO"
    magNumber=None
    #magKey="MAG_APER"
    #magErrKey="MAGERR_APER"
    #magNumber=2
        
    # Dust correction setup
    corrDict={'u': 5.155, 'g': 3.793, 'r': 2.751, 'i': 2.086, 'z': 1.479, 'Ks': 0.367}  
    EBMinusVList=[]
    if fieldIDKey == None:
        EBMinusV=getEBMinusV(np.mean(tab['RADeg']), np.mean(tab['decDeg']), optionsDict = optionsDict)
        EBMinusVList.append({'RADeg': np.mean(tab['RADeg']), 'decDeg': np.mean(tab['decDeg']), 'EBMinusV': EBMinusV})
    else:
        fieldNames=np.unique(tab['field'])
        for f in fieldNames:
            mask=np.where(tab['field'] == f)
            EBMinusV=getEBMinusV(np.mean(tab['RADeg'][mask]), np.mean(tab['decDeg'][mask]), optionsDict = optionsDict)
            EBMinusVList.append({'RADeg': np.mean(tab['RADeg'][mask]), 'decDeg': np.mean(tab['decDeg'][mask]), 'EBMinusV': EBMinusV})

    # Work out available bands
    acceptableBands=['u', 'g', 'r', 'i', 'z', 'Ks']
    tabBands=[]
    for key in tab.keys():
        bits=key.split("_")
        if len(bits) > 0 and bits[0] in acceptableBands and bits[0] not in tabBands:
            tabBands.append(bits[0])
   
    # Parse into zCluster format, apply dust correction too
    catalog=[]
    for row in tab: 
        photDict={}
        photDict['id']=row['ID']
        photDict['RADeg']=row['RADeg']
        photDict['decDeg']=row['decDeg']
        rMin=1e6
        for EBMinusVDict in EBMinusVList:
            rDeg=astCoords.calcAngSepDeg(photDict['RADeg'], photDict['decDeg'], EBMinusVDict['RADeg'], EBMinusVDict['decDeg'])
            if rDeg < rMin:
                rMin=rDeg
                EBMinusV=EBMinusVDict['EBMinusV']
        for b in tabBands:
            dustCorrMag=EBMinusV*corrDict[b]
            if magNumber == None:
                photDict[b]=row['%s_%s' % (b, magKey)]-dustCorrMag
                photDict[b+"Err"]=row['%s_%s' % (b, magErrKey)]
            else:
                photDict[b]=row['%s_%s' % (b, magKey)][magNumber]-dustCorrMag
                photDict[b+"Err"]=row['%s_%s' % (b, magErrKey)][magNumber]
        catalog.append(photDict)

    # Get SDSS catalog, if we want to include missing bands
    if 'addSDSS' in optionsDict.keys() and optionsDict['addSDSS'] == True:
        SDSSCat, SDSSFilterCodes=SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = halfBoxSizeDeg, DR = 12)
        RAs_SDSS=[]
        decs_SDSS=[]
        for objDict in SDSSCat:
            RAs_SDSS.append(objDict['RADeg'])
            decs_SDSS.append(objDict['decDeg'])
        RAs_SDSS=np.array(RAs_SDSS)
        decs_SDSS=np.array(decs_SDSS)
        crossMatchRadiusDeg=1.0/3600.0
        SDSSMatched=np.zeros(len(catalog))  # For sanity checking...
        count=0
        for objDict in catalog:
            rDeg=astCoords.calcAngSepDeg(objDict['RADeg'], objDict['decDeg'], RAs_SDSS, decs_SDSS)
            if rDeg.min() < crossMatchRadiusDeg:
                i=np.where(rDeg == rDeg.min())[0][0]
                SDSSObjDict=SDSSCat[i]
                keysToAdd=['u', 'g', 'r', 'i', 'z']
                for k in keysToAdd:
                    if k not in objDict.keys():
                        objDict['%s' % (k)]=SDSSObjDict[k]
                        objDict['%sErr' % (k)]=SDSSObjDict["%sErr" % (k)]
                    # For sanity check below - we will do this properly somewhere else
                    #objDict['SDSS_%s' % (k)]=SDSSObjDict[k]
                    #objDict['SDSS_%sErr' % (k)]=SDSSObjDict["%sErr" % (k)]
                SDSSMatched[count]=1
            count=count+1
        # Sanity checking
        #mag_cat=[]
        #mag_sdss=[]
        #for i in SDSSMatched.nonzero()[0]:
            #mag_cat.append(catalog[i]['z'])
            #mag_sdss.append(catalog[i]['SDSS_z'])
        #mag_cat=np.array(mag_cat)
        #mag_sdss=np.array(mag_sdss)
        #mask=np.less(mag_sdss, 18)
        #print "SDSS vs FITS table galaxy mag sanity check"
        #IPython.embed()
        #sys.exit()
        
    # This guards against us letting through an empty catalog -  e.g. if not in Stripe 82, we get a 
    # different run number, so object doesn't get added to catalog.
    if catalog == []:
        catalog=None
    
    return catalog

#-------------------------------------------------------------------------------------------------------------
def FITSRetriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, optionsDict = {}):
    """Parses a FITS catalog made by e.g., soi_makecatalogs.py.
    
    NOTE: No star-galaxy separation is applied in this currently.
    
    Here optionsDict needs to include 'fileName' key.
    
    If 'addSDSS': True in optionsDict, then we fetch an SDSS catalog at the location, cross match it against 
    the FITS catalog, and add in info for bands we can't find (e.g., g-band). We don't add SDSS objects which
    aren't detected in the FITS catalog.

    """
    
    try:
        tab=atpy.Table().read(optionsDict['fileName'])
    except:
        raise Exception, "assumed database is a FITS table file, but failed to read"
    
    # If the position isn't actually in our galaxy catalog, we give up now
    rDeg=astCoords.calcAngSepDeg(RADeg, decDeg, tab['RADeg'], tab['decDeg'])
    if rDeg.min() > halfBoxSizeDeg:
        print "... no galaxies found in FITS catalog near RA, dec = (%.6f, %.6f) ..." % (RADeg, decDeg)
        return None
    tab=tab[np.where(rDeg < halfBoxSizeDeg)]
    
    catalog=parseFITSPhotoTable(tab, optionsDict = optionsDict)

    return catalog
