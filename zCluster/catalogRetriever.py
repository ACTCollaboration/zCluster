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
import astropy.table as atpy
import astropy.io.fits as pyfits
import astropy.io.ascii as pyascii
import re
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
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
        os.makedirs(CACHE_DIR, exist_ok = True)

#-------------------------------------------------------------------------------------------------------------
def addWISEPhotometry(RADeg, decDeg, catalog, halfBoxSizeDeg = 36.0/60.):
    """This is an option that can be enabled in other retriever functions by adding 'addWISE': True in
    optionsDict. For this to work, the passbandSet used by PhotoRedshiftEngine must also have the WISE bands
    defined. We can probably tidy this up a bit...
       
    """
    
    cacheDir="WISECache"
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)
    
    # Note that using size= in query below is broken (even though it is on IPAC docs)
    outFileName=cacheDir+os.path.sep+"unWISE_%.6f_%.6f.vot" % (RADeg, decDeg)
    gotFileSuccessfully=False
    while gotFileSuccessfully == False:
        if os.path.exists(outFileName) == False:        
            rah, ram, ras=astCoords.decimal2hms(RADeg, ":").split(":")
            dd, dm, ds=astCoords.decimal2dms(decDeg, ":").split(":")
            radiusSize=halfBoxSizeDeg*3600.0#np.degrees(4./astCalc.da(1.0))*3600.0
            # NOTE: URL modified March 2017 as IPAC seems to have changed name from wise_allwise_p3as_psd to allwise_p3as_psd
            url="http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=cone&catalog=allwise_p3as_psd&objstr=%sh+%sm+%ss+%sd+%sm+%ss&radius=%d&outfmt=3" % (rah, ram, ras, dd, dm, ds, int(round(radiusSize))) 
            urllib.request.urlretrieve(url, outFileName)
        
        try:
            tab=atpy.Table().read(outFileName)
            gotFileSuccessfully=True
        except:
            print("... error reading WISE catalog: retrying download in 1 min ...")
            ##return None
            #IPython.embed()
            #sys.exit()
            os.remove(outFileName)
            time.sleep(60)
    
    # Now need to merge tab with catalog
    matchRadiusDeg=2.0/3600.0
    bands=[]
    for key in catalog[0].keys():
        if key not in ['RADeg', 'decDeg', 'id'] and key.find("Err") == -1:
            bands.append(key)
    catRAs=[]
    catDecs=[]
    for objDict in catalog:
        catRAs.append(objDict['RADeg'])
        catDecs.append(objDict['decDeg'])
        # Defaults: assume not detected
        objDict['w1']=99.0
        objDict['w1Err']=99.0
        objDict['w2']=99.0
        objDict['w2Err']=99.0
    catRAs=np.array(catRAs)
    catDecs=np.array(catDecs)
    objsToAdd=[]
    count=len(catalog)
    for row in tab:
        rDeg=astCoords.calcAngSepDeg(row['ra'], row['dec'], catRAs, catDecs)
        if rDeg.min() < matchRadiusDeg:
            objDict=catalog[np.argmin(rDeg)]
            objDict['w1']=row['w1mpro']+2.699   # Vega -> AB http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
            objDict['w1Err']=row['w1sigmpro']
            objDict['w2']=row['w2mpro']+3.339
            objDict['w2Err']=row['w2sigmpro']
        else:
            count=count+1
            objDict={}
            objDict['id']=count
            objDict['RADeg']=row['ra']
            objDict['decDeg']=row['dec']
            objDict['w1']=row['w1mpro']+2.699
            objDict['w1Err']=row['w1sigmpro']
            objDict['w2']=row['w2mpro']+3.339
            objDict['w2Err']=row['w2sigmpro']
            for b in bands:
                objDict[b]=99.0
                objDict['%sErr' % (b)]=99.0
            objsToAdd.append(objDict)
    catalog=catalog+objsToAdd

    return catalog

#-------------------------------------------------------------------------------------------------------------
def S82Retriever(RADeg, decDeg, halfBoxSizeDeg = 20.2/60.0, optionsDict = {}):
    """Retrieves SDSS Stripe 82 photometry at the given position.
    
    """

    makeCacheDir()
    
    # We used PhotoObj for slit masks, because Galaxy misses some stuff, particularly at high-z
    #tableName="PhotoObj"
    tableName="Galaxy"
    
    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)
        
    url = 'http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp'
            
    outFileName=cacheDir+os.path.sep+"S82_%.4f_%.4f_%.4f_%s.csv" % (RADeg, decDeg, halfBoxSizeDeg, tableName)
    print("... getting SDSS Stripe82 photometry (file: %s) ..." % (outFileName))
            
    if os.path.exists(outFileName) == False:
    
        print("... fetching from the internet ...")

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
        
        # Old, conservative query that ends up missing a lot of galaxies that we can see in S82 i-band
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
        
        # Bail if no chance this is anywhere near S82
        if decMin > 2.5 or decMax < -2.5:
            print("... not close to S82 region ...")
            return None
        
        # Filter SQL so that it'll work
        fsql = ''
        for line in sql.split('\n'):
            fsql += line.split('--')[0] + ' ' + os.linesep;
        
        params=urllib.parse.urlencode({'cmd': fsql, 'format': "csv"})
        response=None
        while response == None:
            try:
                response=urllib.request.urlopen(url+'?%s' % (params))
            except:
                print("Network down? Waiting 30 sec...")
                time.sleep(30)
            
        lines=response.read()
        lines=lines.decode()
        lines=lines.split("\n")

        outFile=open(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=open(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    
    # Parse .csv into catalog
    if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
        catalog=None
    elif lines[1] == 'ERROR\n':
        raise Exception("Error returned by SDSS query - probably too many objects?")
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
                        print("... exceeded server queries per minute limit - waiting ...")
                        time.sleep(70)
                        os.remove(outFileName)
                        return "retry"
                    else:
                        print("what?")
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
                if 'maxMagError' in list(optionsDict.keys()):
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
def DESY3Retriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves DES Y3 photometry at the given position. This assumes you have easyaccess installed
    (https://pypi.python.org/pypi/easyaccess/1.0.7) and access rights for the non-public tables
    in the DES Oracle database.
    
    """
    
    return DESRetriever(RADeg, decDeg, DR = 'Y3', halfBoxSizeDeg = halfBoxSizeDeg, optionsDict = optionsDict)

#-------------------------------------------------------------------------------------------------------------
def DESY3WISERetriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves DES Y3 photometry joined to AllWISE at the given position. This assumes you have easyaccess 
    installed (https://pypi.python.org/pypi/easyaccess/1.0.7) and access rights for the non-public tables
    in the DES Oracle database.
    
    """
    
    return DESRetriever(RADeg, decDeg, DR = 'Y3+WISE', halfBoxSizeDeg = halfBoxSizeDeg, optionsDict = optionsDict)

#-------------------------------------------------------------------------------------------------------------
def DESDR1Retriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves DES DR1 photometry at the given position. This assumes you have easyaccess installed
    (https://pypi.python.org/pypi/easyaccess/1.0.7) and have registered for a login to access the DES
    Oracle database
    
    """
    
    return DESRetriever(RADeg, decDeg, DR = 'DR1', halfBoxSizeDeg = halfBoxSizeDeg, optionsDict = optionsDict)

#-------------------------------------------------------------------------------------------------------------
def DESRetriever(RADeg, decDeg, DR = 'DR1', halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves DES photometry. Assumes you have easyaccess installed
    (https://pypi.python.org/pypi/easyaccess/1.0.7) and the necessary login to access either public or
    proprietary tables in the DES Oracle database.
    
    Use DR = 'DR1' to access the public photometry, DR = 'Y3' to access the current 'Gold' photometry
    
    """

    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)
    
    try:
        connection=optionsDict['connection']
    except:
        raise Exception("No DES database connection")
    
    if decDeg > 5:
        print("... outside DES area - skipping ...")
        return None
    
    outFileName=cacheDir+os.path.sep+"DES%s_%.4f_%.4f_%.4f.csv" % (DR, RADeg, decDeg, halfBoxSizeDeg)      

    # This bit has to be outside, in case we load data from cache
    if DR == 'Y3' or DR == 'Y3+WISE':
        magKey="SOF_CM_MAG_CORRECTED_$BAND"
        magErrKey="SOF_CM_MAG_ERR_$BAND"
    elif DR == 'DR1':
        magKey="MAG_AUTO_$BAND_DERED"
        magErrKey="MAGERR_AUTO_$BAND"
    else:
        raise Exception("didn't recognise requested DES DR")

    if os.path.exists(outFileName) == False or 'refetch' in list(optionsDict.keys()) and optionsDict['refetch'] == True:
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        # Y3 Gold v1.0
        #query="SELECT coadd_object_id, ra, dec, ngmix_cm_mag_g - (3.186 * EBV_SFD98) AS cm_mag_g, ngmix_cm_mag_r - (2.140 * EBV_SFD98) AS cm_mag_r, ngmix_cm_mag_i - (1.569 * EBV_SFD98) AS cm_mag_i, ngmix_cm_mag_z - (1.196 * EBV_SFD98) AS cm_mag_z, ngmix_cm_mag_err_g AS cm_mag_err_g, ngmix_cm_mag_err_r AS cm_mag_err_r, ngmix_cm_mag_err_i AS cm_mag_err_i, ngmix_cm_mag_err_z AS cm_mag_err_z FROM y3_gold_1_0 SAMPLE WHERE ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f AND flag_gold = 0 AND flag_footprint = 1 AND flag_foreground = 0 AND extended_class_mash BETWEEN 3 AND 4" % (RAMin, RAMax, decMin, decMax) 
        # Y3 Gold v2.2
        if DR == 'Y3':
            query="SELECT COADD_OBJECT_ID, RA, DEC, DNF_ZMC_SOF, BPZ_ZMC_SOF, SOF_CM_MAG_CORRECTED_G, SOF_CM_MAG_CORRECTED_R, SOF_CM_MAG_CORRECTED_I, SOF_CM_MAG_CORRECTED_Z, SOF_CM_MAG_ERR_G, SOF_CM_MAG_ERR_R, SOF_CM_MAG_ERR_I, SOF_CM_MAG_ERR_Z FROM Y3_GOLD_2_2 WHERE FLAGS_FOOTPRINT = 1 and FLAGS_FOREGROUND = 0 and bitand(FLAGS_GOLD, 62) = 0 and EXTENDED_CLASS_MASH_SOF = 3 and SOF_CM_MAG_I between 16 and 24 AND RA BETWEEN %.6f AND %.6f AND DEC BETWEEN %.6f and %.6f" % (RAMin, RAMax, decMin, decMax)
        elif DR == 'DR1':
            query="SELECT RA, DEC, MAG_AUTO_G_DERED, MAG_AUTO_R_DERED, MAG_AUTO_I_DERED, MAG_AUTO_Z_DERED, MAGERR_AUTO_G, MAGERR_AUTO_R, MAGERR_AUTO_I, MAGERR_AUTO_Z FROM DR1_MAIN WHERE WAVG_SPREAD_MODEL_I + 3.0*WAVG_SPREADERR_MODEL_I > 0.005 and WAVG_SPREAD_MODEL_I + 1.0*WAVG_SPREADERR_MODEL_I > 0.003 and WAVG_SPREAD_MODEL_I - 1.0*WAVG_SPREADERR_MODEL_I > 0.001 and WAVG_SPREAD_MODEL_I > -1 and IMAFLAGS_ISO_I = 0 and MAG_AUTO_I < 24 and RA BETWEEN %.6f AND %.6f AND DEC BETWEEN %.6f and %.6f" % (RAMin, RAMax, decMin, decMax)
        elif DR == 'Y3+WISE':
            query="SELECT D.COADD_OBJECT_ID, D.RA, D.DEC, DNF_ZMC_SOF, BPZ_ZMC_SOF, SOF_CM_MAG_CORRECTED_G, SOF_CM_MAG_CORRECTED_R, SOF_CM_MAG_CORRECTED_I, SOF_CM_MAG_CORRECTED_Z, SOF_CM_MAG_ERR_G, SOF_CM_MAG_ERR_R, SOF_CM_MAG_ERR_I, SOF_CM_MAG_ERR_Z, W1MPRO, W1SIGMPRO, W2MPRO, W2SIGMPRO FROM Y3_GOLD_2_2 D LEFT OUTER JOIN DES_ADMIN.Y3A2_WISE_DES W ON D.COADD_OBJECT_ID = W.COADD_OBJECT_ID WHERE FLAGS_FOOTPRINT = 1 and FLAGS_FOREGROUND = 0 and bitand(FLAGS_GOLD, 62) = 0 and EXTENDED_CLASS_MASH_SOF = 3 and SOF_CM_MAG_I between 16 and 24 AND D.RA BETWEEN %.6f AND %.6f AND D.DEC BETWEEN %.6f and %.6f" % (RAMin, RAMax, decMin, decMax)
        else:
            raise Exception("didn't recognise requested DES DR")

        # This should restart the connection if it drops
        if connection.ping() == False:
            print("... DES database connection lost: reconnecting ...")
            import easyaccess as ea
            connection=ea.connect()
        if connection.ping() == True:
            connection.query_and_save(query, outFileName)
        else:
            raise Exception("DES database connection lost")

    # No output is written if the query returns nothing... so we'll write a blank file to save repeating the query next time
    catalog=[]
    if os.path.exists(outFileName) == False:
        outFile=open(outFileName, "w")
        outFile.write("# No objects returned by query\n")
        outFile.close()
    else:
        
        tab=atpy.Table().read(outFileName, format = 'csv')

        idCount=0
        for row in tab:
            idCount=idCount+1
            photDict={}
            photDict['id']=idCount    # just so we have something - we could use COADD_OBJECT_ID but skipping for now
            photDict['RADeg']=row['RA']
            photDict['decDeg']=row['DEC']
            photDict['g']=row[magKey.replace("$BAND", "G")]
            photDict['r']=row[magKey.replace("$BAND", "R")]
            photDict['i']=row[magKey.replace("$BAND", "I")]
            photDict['z']=row[magKey.replace("$BAND", "Z")]
            photDict['gErr']=row[magErrKey.replace("$BAND", "G")]
            photDict['rErr']=row[magErrKey.replace("$BAND", "R")]
            photDict['iErr']=row[magErrKey.replace("$BAND", "I")]
            photDict['zErr']=row[magErrKey.replace("$BAND", "Z")]
            if DR.find("+WISE") != -1:
                # Convert to AB mags and handle missing values
                if row['W1MPRO'] < 0:
                    w1Mag=99.0
                    w1MagErr=99.0
                else:
                    w1Mag=row['W1MPRO']+2.699
                    w1MagErr=row['W1SIGMPRO']
                if row['W2MPRO'] < 0:
                    w2Mag=99.0
                    w2MagErr=99.0
                else:
                    w2Mag=row['W2MPRO']+3.339
                    w2MagErr=row['W2SIGMPRO']
                photDict['w1']=w1Mag
                photDict['w1Err']=w1MagErr
                photDict['w2']=w2Mag
                photDict['w2Err']=w2MagErr
            
            # Apply mag error cuts if given
            # For PS1, missing values are -999 - our current checkMagErrors routine will fish those out
            # We're just making the mag unconstrained here (missing data), rather than applying a limit
            # If we don't have a minimum of three useful bands, reject
            if 'maxMagError' in list(optionsDict.keys()):
                keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['g', 'r', 'i', 'z'])
            else:
                keep=True
            
            ## Additional colour cuts (as used in some DES papers - see e.g., splashback draft)
            #gr=photDict['g']-photDict['r']
            #ri=photDict['r']-photDict['i']
            #iz=photDict['i']-photDict['z']
            #if keep == True and -1 < gr < 3 and -1 < ri < 2.5 and -1 < iz < 2:
                #keep=True
            #else:
                #keep=False
                
            if keep == True:
                catalog.append(photDict)
    
    #if 'addWISE' in optionsDict.keys() and optionsDict['addWISE'] == True:
        #catalog=addWISEPhotometry(RADeg, decDeg, catalog, halfBoxSizeDeg = halfBoxSizeDeg)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def KiDSDR4Retriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves KiDS DR4 photometry (which includes VIKING IR bands).
    
    """
    
    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)
        
    # This defines the rough final KIDS area (split into two fields)
    inKIDSRegion=False
    if RADeg > 120 and RADeg < 240 and decDeg > -5 and decDeg < 5:
        inKIDSRegion=True
    if RADeg > 330 and decDeg > -37 and decDeg < -25:
        inKIDSRegion=True
    if RADeg < 60 and decDeg > -37 and decDeg < -25:
        inKIDSRegion=True
    if inKIDSRegion == False:
        print("... outside KIDS area - skipping ...")
        return None

    outFileName=cacheDir+os.path.sep+"KiDSDR4_%.4f_%.4f_%.4f.fits" % (RADeg, decDeg, halfBoxSizeDeg)      
    if os.path.exists(outFileName) == False or 'refetch' in list(optionsDict.keys()) and optionsDict['refetch'] == True:
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        print("... downloading catalog %s ..." % (outFileName))
        query="select id, RAJ2000, DECJ2000, MAG_GAAP_u, MAG_GAAP_g, MAG_GAAP_r, MAG_GAAP_i, MAG_GAAP_Z, MAG_GAAP_Y, MAG_GAAP_J, MAG_GAAP_H, MAG_GAAP_Ks, MAGERR_GAAP_u, MAGERR_GAAP_g, MAGERR_GAAP_r, MAGERR_GAAP_i, MAGERR_GAAP_Z, MAGERR_GAAP_Y, MAGERR_GAAP_J, MAGERR_GAAP_H, MAGERR_GAAP_Ks, SG2DPHOT from KiDS_DR4_0_ugriZYJHKs_cat_fits_V3 where RAJ2000 BETWEEN %.6f AND %.6f AND DECJ2000 BETWEEN %.6f and %.6f and SG2DPHOT = 0;" % (RAMin, RAMax, decMin, decMax)        
        r=requests.get('http://archive.eso.org/tap_cat/sync?', 
                        params = {'REQUEST': 'doQuery',
                                  'LANG': 'ADQL',
                                  'MAXREC': 1000000,
                                  'FORMAT': 'fits',
                                  'QUERY': query},
                        stream = True)
        with open(outFileName, 'wb') as outFile: 
            r.raw.decode_content = True
            outFile.write(r.raw.data) 
    
    # Load/parse table
    print("... reading catalog %s ..." % (outFileName))
    try:
        tab=atpy.Table().read(outFileName)
    except ValueError:
        print("... no objects in catalog - skipping ...")
        return None
    
    magKey="MAG_GAAP_$BAND"
    magErrKey="MAGERR_GAAP_$BAND"
    
    # First, get rid of nans or nonsensical values
    bands=['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']
    for b in bands:
        #tab=tab[np.where(np.isnan(tab[magKey.replace("$BAND", b)]) == False)]
        #tab=tab[np.where(np.isnan(tab[magErrKey.replace("$BAND", b)]) == False)]
        tab=tab[np.where(tab[magKey.replace("$BAND", b)] > 0)]
        tab=tab[np.where(tab[magErrKey.replace("$BAND", b)] > 0)]
        
    # NOTE: Need to add JHKs mags
    idCount=0
    catalog=[]
    for row in tab:
        idCount=idCount+1
        photDict={}
        photDict['id']=idCount    # just so we have something - we could use ID but skipping for now
        photDict['RADeg']=row['RAJ2000']
        photDict['decDeg']=row['DECJ2000']
        # KiDS DR4 mags we're using are already extinction corrected
        photDict['u']=row[magKey.replace("$BAND", "u")]#-row['EXT_SFD_U']
        photDict['g']=row[magKey.replace("$BAND", "g")]#-row['EXT_SFD_G']
        photDict['r']=row[magKey.replace("$BAND", "r")]#-row['EXT_SFD_R']
        photDict['i']=row[magKey.replace("$BAND", "i")]#-row['EXT_SFD_I']
        photDict['Z']=row[magKey.replace("$BAND", "Z")]#-row['EXT_SFD_U']
        photDict['Y']=row[magKey.replace("$BAND", "Y")]#-row['EXT_SFD_G']
        photDict['J']=row[magKey.replace("$BAND", "J")]#-row['EXT_SFD_R']
        photDict['H']=row[magKey.replace("$BAND", "H")]#-row['EXT_SFD_I']
        photDict['Ks']=row[magKey.replace("$BAND", "Ks")]#-row['EXT_SFD_I']
        photDict['uErr']=row[magErrKey.replace("$BAND", "u")]
        photDict['gErr']=row[magErrKey.replace("$BAND", "g")]
        photDict['rErr']=row[magErrKey.replace("$BAND", "r")]
        photDict['iErr']=row[magErrKey.replace("$BAND", "i")]
        photDict['ZErr']=row[magErrKey.replace("$BAND", "Z")]
        photDict['YErr']=row[magErrKey.replace("$BAND", "Y")]
        photDict['JErr']=row[magErrKey.replace("$BAND", "J")]
        photDict['HErr']=row[magErrKey.replace("$BAND", "H")]
        photDict['KsErr']=row[magErrKey.replace("$BAND", "Ks")]
        if 'maxMagError' in list(optionsDict.keys()):
            keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['u', 'g', 'r', 'i'])
        else:
            keep=True
                    
        if keep == True:
            catalog.append(photDict)
        
    return catalog
    
#-------------------------------------------------------------------------------------------------------------
def ATLASDR4Retriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves VST ATLAS DR4 photometry via ESO.
    
    """

    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)
        
    outFileName=cacheDir+os.path.sep+"ATLASDR4_%.4f_%.4f_%.4f.fits" % (RADeg, decDeg, halfBoxSizeDeg)          
    if os.path.exists(outFileName) == False or 'refetch' in list(optionsDict.keys()) and optionsDict['refetch'] == True:
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        print("... downloading catalog %s ..." % (outFileName)) 
        query="select sourceID, ra2000, dec2000, uPetroMag, uPetroMagErr, gPetroMag, gPetroMagErr, rPetroMag, rPetroMagErr, iPetroMag, iPetroMagErr, zPetroMag, zPetroMagErr, uAperMagNoAperCorr3, uAperMag3Err, gAperMagNoAperCorr3, gAperMag3Err, rAperMagNoAperCorr3, rAperMag3Err, iAperMagNoAperCorr3, iAperMag3Err, zAperMagNoAperCorr3, zAperMag3Err, aU, aG, aR, aI, aZ, pGalaxy, pStar, pNoise, pSaturated from atlas_er4_ugriz_catMetaData_fits_V3 where ra2000 BETWEEN %.4f AND %.4f AND dec2000 BETWEEN %.4f and %.4f and priOrSec = 0;" % (RAMin, RAMax, decMin, decMax)   
        r=requests.get('http://archive.eso.org/tap_cat/sync?', 
                        params = {'REQUEST': 'doQuery',
                                  'LANG': 'ADQL',
                                  'MAXREC': 1000000,
                                  'FORMAT': 'fits',
                                  'QUERY': query},
                        stream = True)
        with open(outFileName, 'wb') as outFile: 
            r.raw.decode_content = True
            outFile.write(r.raw.data) 
            
    # Load/parse table
    print("... reading catalog %s ..." % (outFileName))
    tab=atpy.Table().read(outFileName)
    if len(tab) == 0:
        print("... no objects in catalog - skipping ...")
        return None

    magKey="$BANDPETROMAG"
    magErrKey="$BANDPETROMAGERR"
    #magKey="$BANDAPERMAGNOAPERCORR3"
    #magErrKey="$BANDAPERMAG3ERR"
    
    # First, get rid of nans or nonsensical values
    # We apply extinction correction in-place here
    # NOTE: Zapping u (most uncertain with regards photo calib I think) makes no significant difference
    bands=['u', 'g', 'r', 'i', 'z']
    for b in bands:
        tab[magKey.replace("$BAND", b.upper())][np.isnan(tab[magKey.replace("$BAND", b.upper())])]=99.0
        tab[magErrKey.replace("$BAND", b.upper())][np.isnan(tab[magErrKey.replace("$BAND", b.upper())])]=99.0
        tab[magKey.replace("$BAND", b.upper())]=tab[magKey.replace("$BAND", b.upper())]-tab['A%s' % (b.upper())]
        
    # Transform to SDSS in-place - maybe this works better?
    # These transforms are from doing some algebra on the ones in the VST-ATLAS paper and throwing away ~0.00x terms
    # NOTE: These don't work at all well (much worse than not doing this)
    #iSDSS=tab[magKey.replace("$BAND", 'I')]-0.025
    #zSDSS=(tab[magKey.replace("$BAND", 'Z')]-0.04*tab[magKey.replace("$BAND", 'I')]+0.04)/0.96
    #rSDSS=(tab[magKey.replace("$BAND", 'R')]+0.0316*tab[magKey.replace("$BAND", 'G')])/1.03
    #gSDSS=(tab[magKey.replace("$BAND", 'G')]/0.95)-0.048*(tab[magKey.replace("$BAND", 'R')]+0.316*tab[magKey.replace("$BAND", 'G')])-0.063
    #uSDSS=1.01*(tab[magKey.replace("$BAND", 'U')]-0.01*gSDSS+0.27)
    #tab[magKey.replace("$BAND", 'U')]=uSDSS
    #tab[magKey.replace("$BAND", 'G')]=gSDSS
    #tab[magKey.replace("$BAND", 'R')]=rSDSS
    #tab[magKey.replace("$BAND", 'I')]=iSDSS
    #tab[magKey.replace("$BAND", 'Z')]=zSDSS
    
    # Empirical corrections (based on own comparison with SDSS galaxy photometry for small subsample)
    # Again, these don't help at all really
    #tab[magKey.replace("$BAND", 'U')]=tab[magKey.replace("$BAND", 'U')]-0.286
    #tab[magKey.replace("$BAND", 'G')]=tab[magKey.replace("$BAND", 'G')]-0.080
    #tab[magKey.replace("$BAND", 'R')]=tab[magKey.replace("$BAND", 'R')]-0.009
    #tab[magKey.replace("$BAND", 'I')]=tab[magKey.replace("$BAND", 'I')]-0.084
    #tab[magKey.replace("$BAND", 'Z')]=tab[magKey.replace("$BAND", 'Z')]-0.073
    
    # Classification cuts
    tab=tab[np.where(tab['PGALAXY'] > 0.5)]
    tab=tab[np.where(tab['PNOISE'] < 0.01)]
    tab=tab[np.where(tab['PSATURATED'] < 0.01)]
    
    idCount=0
    catalog=[]
    for row in tab:
        idCount=idCount+1
        photDict={}
        photDict['id']=idCount    # just so we have something - we could use ID but skipping for now
        photDict['RADeg']=row['RA2000']
        photDict['decDeg']=row['DEC2000']
        for b in bands:
            photDict[b]=row[magKey.replace("$BAND", b.upper())]
            photDict[b+'Err']=row[magErrKey.replace("$BAND", b.upper())]
        if 'maxMagError' in list(optionsDict.keys()):
            keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = bands)
        else:
            keep=True
        if keep == True:
            catalog.append(photDict)
        
    return catalog

#------------------------------------------------------------------------------------------------------------
def fixcolnames(tab):
    """Fix column names returned by the casjobs query. Needed by PS1Retriever.
    
    Parameters
    ----------
    tab (astropy.table.Table): Input table

    Returns reference to original table with column names modified"""

    pat = re.compile(r'\[(?P<name>[^[]+)\]')
    for c in tab.colnames:
        m = pat.match(c)
        if not m:
            raise ValueError("Unable to parse column name '{}'".format(c))
        newname = m.group('name')
        tab.rename_column(c,newname)
    return tab

#-------------------------------------------------------------------------------------------------------------
def PS1Retriever(RADeg, decDeg, halfBoxSizeDeg = 25.5/60.0, optionsDict = {}):
    """Retrieves PS1 photometry at the given position.
        
    """
    
    import mastcasjobs

    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)

    outFileName=cacheDir+os.path.sep+"PS1_%.4f_%.4f_%.4f.fits" % (RADeg, decDeg, halfBoxSizeDeg)      
    print("... getting PS1 photometry (file: %s) ..." % (outFileName))

    if decDeg < -30:
        print("... outside PS1 area - skipping ...")
        return None
    
    if os.path.exists(outFileName) == False or 'refetch' in list(optionsDict.keys()) and optionsDict['refetch'] == True:
        print("... fetching from the internet ...")
        query="""select o.objID, o.raMean, o.decMean,
        m.gMeanKronMag, m.rMeanKronMag, m.iMeanKronMag, m.zMeanKronMag, m.yMeanKronMag, m.gMeanKronMagErr, m.rMeanKronMagErr, m.iMeanKronMagErr, m.zMeanKronMagErr, m.yMeanKronMagErr
        from fGetNearbyObjEq(%.6f, %.6f, %.6f) nb
        inner join ObjectThin o on o.objid=nb.objid and o.nDetections > 3
        inner join MeanObject m on o.objid=m.objid and o.uniquePspsOBid=m.uniquePspsOBid""" % (RADeg, decDeg, halfBoxSizeDeg*60)
        jobs=optionsDict['jobs']
        #jobs=mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
        results=jobs.quick(query, task_name="python cone search")
        tab=fixcolnames(pyascii.read(results))
        tab.write(outFileName, overwrite = True)
    else:
        tab=atpy.Table().read(outFileName)
        
    # Parse table into catalog
    if len(tab) == 0:
        catalog=None
    else:
        EBMinusV=getEBMinusV(RADeg, decDeg, optionsDict = optionsDict) # assume same across field
        catalog=[]
        idCount=0
        bands=['g', 'r', 'i', 'z', 'y']
        for row in tab:
            idCount=idCount+1
            photDict={}
            photDict['id']=idCount    # just so we have something - we could use PS1 ID but skipping for now
            photDict['RADeg']=row['raMean']
            photDict['decDeg']=row['decMean']
            for b in bands:
                photDict[b]=row['%sMeanKronMag' % (b)]
                photDict['%sErr' % (b)]=row['%sMeanKronMagErr' %  (b)]

            # Correct for dust extinction
            # Taken from: http://www.mso.anu.edu.au/~brad/filters.html
            photDict['g']=photDict['g']-EBMinusV*3.322
            photDict['r']=photDict['r']-EBMinusV*2.544 
            photDict['i']=photDict['i']-EBMinusV*2.265 
            photDict['z']=photDict['z']-EBMinusV*1.846 
            photDict['y']=photDict['y']-EBMinusV*1.570 

            # Apply mag error cuts if given
            # For PS1, missing values are -999 - our current checkMagErrors routine will fish those out
            # We're just making the mag unconstrained here (missing data), rather than applying a limit
            # If we don't have a minimum of three useful bands, reject
            if 'maxMagError' in list(optionsDict.keys()):
                keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['g', 'r', 'i', 'z'], minBands = 4)
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
    
    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)
    
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
                                    
    print("... getting SDSS DR%d photometry (file: %s) ..." % (DR, outFileName))
        
    if os.path.exists(outFileName) == False or 'refetch' in list(optionsDict.keys()) and optionsDict['refetch'] == True:
        
        print("... fetching from the internet ...")
        
        # Clean galaxy photometry query - note flags for r-band only, may want to change
        # We may want to add something that does multiple queries if we want a bigger area from the standard
        # SDSS query interface
        # NOTE: disabled the 'clean' photometry flags on 21/09/2015
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        # Sanity check: if RAMin, RAMax have order reversed, SDSS gives us nothing
        #if RAMin > RAMax:
            #print "RAMin, RAMax reversed"
            #IPython.embed()
            #sys.exit()
            #newRAMax=RAMin
            #RAMin=RAMax
            #RAMax=newRAMax
        #if decMin > decMax:
            #print "decMin, decMax reversed"
            #IPython.embed()
            #sys.exit()
            #newDecMax=decMin
            #decMin=decMax
            #decMax=newDecMax
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
    
        params=urllib.parse.urlencode({'cmd': fsql, 'format': "csv"})
        response=None
        while response == None:
            try:
                response=urllib.request.urlopen(url+'%s' % (params))
            except:
                print("Network down? Waiting 30 sec...")
                time.sleep(30)

        # Some faffing about here because of python2 -> python3 
        lines=response.read()
        lines=lines.splitlines()        
        outFile=open(outFileName, "w")
        strLines=[]
        for line in lines:
            strLines.append(line.decode('utf-8'))
            outFile.write(strLines[-1]+"\n")
        outFile.close()
        lines=strLines
    
    else:
        
        inFile=open(outFileName, "r")
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
                    print("... problem with file %s - removing and retrying in 70 sec ..." % (outFileName))
                    os.remove(outFileName)
                    time.sleep(70)
                    return "retry"                    
                    #if lines[1][:46] == '"ERROR: Maximum 60 queries allowed per minute.':
                        #print "... exceeded server queries per minute limit - waiting ..."
                        #time.sleep(70)
                        #os.remove(outFileName)
                        #return "retry"
                    #else:
                        #if line.find("<html>") != -1:
                            #time.sleep(70)
                            #os.remove(outFileName)
                            #return "retry"
                        #else:
                            #print "what?"
                            #IPython.embed()
                            #sys.exit()
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
                if 'maxMagError' in list(optionsDict.keys()):
                    keep=checkMagErrors(photDict, optionsDict['maxMagError'])
                else:
                    keep=True
                if keep == True:
                    catalog.append(photDict)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def DECaLSRetriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves DECaLS DR8 tractor catalogs (if they exist) at the given position. Cuts the catalog to the
    radius specified by halfBoxSizeDeg.

    """

    makeCacheDir()
    
    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir, exist_ok = True)

    # Organised such that after this, have subdir with degrees RA (e.g. 000/ is 0 < RADeg < 1 degree)
    basePath="https://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/south/tractor/"

    outFileName=cacheDir+os.path.sep+"DECaLSDR8_%.4f_%.4f_%.2f.fits" % (RADeg, decDeg,
                                                                            halfBoxSizeDeg)
    
    print("... getting DECaLS photometry ...")

    bricksTab=optionsDict['bricksTab']
    DR8Tab=optionsDict['DR8Tab']

    # Find matching tractor catalogs and download/cache .fits table files
    # Previously this would fail for very small search areas - hence added centre coords search
    # NOTE: We've doubled asked for halfBoxSizeDeg here to get all the matching tiles
    # Otherwise we tend to miss...
    RAMask=np.logical_and(np.greater(RADeg, bricksTab['RA1']), np.less(RADeg, bricksTab['RA2']))
    decMask=np.logical_and(np.greater(decDeg, bricksTab['DEC1']), np.less(decDeg, bricksTab['DEC2']))
    centreMask=np.logical_and(RAMask, decMask)
    RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, 2*halfBoxSizeDeg)
    mask=np.logical_and(np.greater(bricksTab['RA1'], RAMin), np.less(bricksTab['RA2'], RAMax))
    mask=np.logical_and(mask, np.greater(bricksTab['DEC1'], decMin))
    mask=np.logical_and(mask, np.less(bricksTab['DEC2'], decMax))
    mask=np.logical_or(mask, centreMask)
    matchTab=bricksTab[np.where(mask)]
    count=0
    tractorTabs=[]
    try:
        matchTab=atpy.join(matchTab, DR8Tab, keys = 'BRICKNAME')
    except:
        # Not in DECaLS?
        print("... no match between bricks and DR8 tables - %s ..." % (outFileName))
        return None
    for mrow in matchTab:
        subDir="%03d" % np.floor(mrow['RA'])
        url=basePath+subDir
        url=url+os.path.sep+"tractor-"+mrow['BRICKNAME']+".fits"
        os.makedirs(cacheDir+os.path.sep+subDir, exist_ok = True)
        fileName=cacheDir+os.path.sep+subDir+os.path.sep+"tractor-%s.fits" % (mrow['BRICKNAME'])
        if os.path.exists(fileName) == False:
            print("... retrieving tractor catalog from web: %s ..." % (url))
            try:
                urllib.request.urlretrieve(url, filename = fileName)
            except:
                with open("wget_failed.sh", "a") as outFile:
                    outFile.write("wget -nc %s\n" % (url))
                print("... WARNING: failed to fetch from %s" % (url))
        try:
            tractorTab=atpy.Table.read(fileName)
            tractorTabs.append(tractorTab)
        except:
            print("... possibly a 404 error for %s - check if cached file is corrupted ..." % (fileName))
    if 'downloadOnly' in optionsDict.keys() and optionsDict['downloadOnly'] == True:
        return None
    
    # Stitch catalogs together
    if len(tractorTabs) > 0:
        tab=atpy.vstack(tractorTabs)
        
        # Remove stars
        tab=tab[np.where(tab['type'] != 'PSF')]
        
        # Cut to asked for size
        rDeg=astCoords.calcAngSepDeg(RADeg, decDeg, tab['ra'], tab['dec'])
        mask=np.less(rDeg, halfBoxSizeDeg)
        tab=tab[mask]
                
        # WISE fluxes are available...
        bands=['g', 'r', 'z', "w1", "w2"]# , 'Y']

        # Convert nanomaggies to mags and do extinction correction
        bricksInTab=np.unique(tab['brickname'])
        for brickName in bricksInTab:
            brickExtTab=DR8Tab[np.where(DR8Tab['BRICKNAME'] == brickName)]
            brickIndices=np.where(tab['brickname'] == brickName)
            brickMask=np.zeros(len(tab), dtype = bool)
            brickMask[brickIndices]=True
            for b in bands:
                magLabel='%s_mag' % (b)
                magErrLabel='%s_magErr' % (b)
                extKey='ext_%s' % (b)
                colsToAdd=[magLabel, magErrLabel]
                for col in colsToAdd:
                    if col not in list(tab.keys()):
                        tab.add_column(atpy.Column(np.ones(len(tab))*99., col))
                # Dust correction: correct flux upwards using mw_transmission (0-1)
                # Instead of using extinction mags below (seems to make no difference)
                nanoMaggies=np.array(tab['flux_%s' % (b)] / tab['mw_transmission_%s' % (b)])
                validMask=np.greater(nanoMaggies, 0.)
                mask=np.logical_and(brickMask, validMask)
                tab[magLabel][mask]=-2.5*np.log10(nanoMaggies[mask])+22.5
                validFluxErr=np.sqrt(1./tab['flux_ivar_%s' % (b)][mask])
                tab[magErrLabel][mask]=1./(nanoMaggies[mask]/validFluxErr)
                # Old...
                #if extKey in list(brickExtTab.keys()):
                    #dustCorrMag=brickExtTab[extKey][0]
                    #tab[magLabel][mask]=tab[magLabel][mask]-dustCorrMag
        
        catalog=[]
        for row in tab: 
            photDict={}
            photDict['id']=row['objid']
            photDict['RADeg']=row['ra']
            photDict['decDeg']=row['dec']
            for b in bands:
                photDict[b]=row['%s_mag' % (b)]
                photDict[b+"Err"]=row['%s_magErr' % (b)]
            if 'maxMagError' in list(optionsDict.keys()):
                keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = bands)
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
    
    if 'altCacheDir' in list(optionsDict.keys()):
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

    print("... getting CFHTLenS photometry (file: %s) ..." % (outFileName)) 
        
    if os.path.exists(outFileName) == False:
        
        print("... fetching from the internet ...")
            
        response=None
        while response == None:
            try:
                response=urllib.request.urlopen(url)
            except:
                print("Network down, or CADC changed URLs again? Waiting 30 sec...")
                time.sleep(30)
                
        lines=response.read()
        lines=lines.split("\n")

        outFile=open(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=open(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()

    # Parse .csv into catalog
    if len(lines) < 10:
        catalog=None
    else:
        catalog=[]
        idCount=0
        # NOTE: Extinction corrections already applied to CFHTLenS catalogs
        for line in lines[1:]: # first line always heading
            if len(line) > 3:
                photDict={}
                idCount=idCount+1
                bits=line.replace("\n", "").split()                
                photDict['id']=idCount    # just so we have something
                try:
                    photDict['RADeg']=float(bits[1])
                except:
                    print("... problem with file %s - removing and retrying in 70 sec ..." % (outFileName))
                    os.remove(outFileName)
                    time.sleep(70)
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
                                                
                # Apply mag error cuts if given
                # We're just making the mag unconstrained here (missing data), rather than applying a limit
                # If we don't have a minimum of three useful bands, reject
                if 'maxMagError' in list(optionsDict.keys()):
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
def getEBMinusV(RADeg, decDeg, optionsDict = {}):
    """Get E(B-V) due to Galactic dust from Schlegel maps, using IRSA web service
    
    """

    makeCacheDir()
    if 'altCacheDir' in list(optionsDict.keys()):
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
        
    fileName=cacheDir+os.path.sep+"SchlegelIRSA_%.6f_%.6f.xml" % (RADeg, decDeg)
    if os.path.exists(fileName) == False:
        url="http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr="+str(RADeg)+"+"+str(decDeg)+"+equ+J2000"
        urllib.request.urlretrieve(url, filename = fileName)
    
    inFile=open(fileName, "r")
    lines=inFile.readlines()
    inFile.close()
    for i in range(len(lines)):
        if lines[i].find("refPixelValueSFD") != -1:
            break
    try:
        EBMinusV=float(lines[i+1].split(" (")[0])
    except:
        raise Exception("failed to parse IRSA E(B-V) value - problem with file %s?" % (fileName))

    return EBMinusV

#-------------------------------------------------------------------------------------------------------------
def parseFITSPhotoTable(tab, fieldIDKey = None, optionsDict = {}):
    """Parse .fits table into catalog list of dictionaries format. Extinction correction due to Galactic
    dust will be applied, using IRSA web service. If fieldIDKey == None (the default), this is done
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
    #corrDict={'u': 5.155, 'g': 3.793, 'r': 2.751, 'i': 2.086, 'z': 1.479, 'Ks': 0.367}  
    #corrDict={'u': 5.155, 'g': 3.793, 'r': 2.751, 'i': 2.086, 'z': 1.479, 'Ks': 0.367}  
    #EBMinusVList=[]
    #if fieldIDKey == None:
        #EBMinusV=getEBMinusV(np.mean(tab['RADeg']), np.mean(tab['decDeg']), optionsDict = optionsDict)
        #EBMinusVList.append({'RADeg': np.mean(tab['RADeg']), 'decDeg': np.mean(tab['decDeg']), 'EBMinusV': EBMinusV})
    #else:
        #fieldNames=np.unique(tab['field'])
        #for f in fieldNames:
            #mask=np.where(tab['field'] == f)
            #EBMinusV=getEBMinusV(np.mean(tab['RADeg'][mask]), np.mean(tab['decDeg'][mask]), optionsDict = optionsDict)
            #EBMinusVList.append({'RADeg': np.mean(tab['RADeg'][mask]), 'decDeg': np.mean(tab['decDeg'][mask]), 'EBMinusV': EBMinusV})

    # Work out available bands
    acceptableBands=['u', 'g', 'r', 'i', 'z', 'J', 'H', 'Ks']
    tabBands=[]
    for key in list(tab.keys()):
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
        #for EBMinusVDict in EBMinusVList:
            #rDeg=astCoords.calcAngSepDeg(photDict['RADeg'], photDict['decDeg'], EBMinusVDict['RADeg'], EBMinusVDict['decDeg'])
            #if rDeg < rMin:
                #rMin=rDeg
                #EBMinusV=EBMinusVDict['EBMinusV']
        #for b in tabBands:
            #dustCorrMag=EBMinusV*corrDict[b]
            #if magNumber == None:
                #photDict[b]=row['%s_%s' % (b, magKey)]-dustCorrMag
                #photDict[b+"Err"]=row['%s_%s' % (b, magErrKey)]
            #else:
                #photDict[b]=row['%s_%s' % (b, magKey)][magNumber]-dustCorrMag
                #photDict[b+"Err"]=row['%s_%s' % (b, magErrKey)][magNumber]
        for b in tabBands:
            if magNumber == None:
                photDict[b]=row['%s_%s' % (b, magKey)]
                photDict[b+"Err"]=row['%s_%s' % (b, magErrKey)]
            else:
                photDict[b]=row['%s_%s' % (b, magKey)][magNumber]
                photDict[b+"Err"]=row['%s_%s' % (b, magErrKey)][magNumber]
        # Optional: spec-zs for fitting for zero point offsets
        if 'z_spec' in tab.keys():
            photDict['z_spec']=row['z_spec']
        if 'maxMagError' in list(optionsDict.keys()):
            keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = acceptableBands)
        else:
            keep=True
        catalog.append(photDict)

    # Get SDSS catalog, if we want to include missing bands
    if 'addSDSS' in list(optionsDict.keys()) and optionsDict['addSDSS'] == True:
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
                    if k not in list(objDict.keys()):
                        objDict['%s' % (k)]=SDSSObjDict[k]
                        objDict['%sErr' % (k)]=SDSSObjDict["%sErr" % (k)]
                    # For sanity check below - we will do this properly somewhere else
                    #objDict['SDSS_%s' % (k)]=SDSSObjDict[k]
                    #objDict['SDSS_%sErr' % (k)]=SDSSObjDict["%sErr" % (k)]
                SDSSMatched[count]=1
            count=count+1
        
    if catalog == []:
        catalog=None

    return catalog

#-------------------------------------------------------------------------------------------------------------
def FITSRetriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
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
        raise Exception("assumed database is a FITS table file, but failed to read")
    
    # If the position isn't actually in our galaxy catalog, we give up now
    rDeg=astCoords.calcAngSepDeg(RADeg, decDeg, tab['RADeg'], tab['decDeg'])
    if rDeg.min() > halfBoxSizeDeg:
        print("... no galaxies found in FITS catalog near RA, dec = (%.6f, %.6f) ..." % (RADeg, decDeg))
        return None
    tab=tab[np.where(rDeg < halfBoxSizeDeg)]
    
    catalog=parseFITSPhotoTable(tab, optionsDict = optionsDict)

    return catalog
