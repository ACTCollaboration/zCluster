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
            print "... not close to S82 region ..."
            return None
        
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
def DESY3Retriever(RADeg, decDeg, halfBoxSizeDeg = 36.0/60.0, optionsDict = {}):
    """Retrieves DES Y3 photometry at the given position. This assumes you have easyaccess installed
    (https://pypi.python.org/pypi/easyaccess/1.0.7) and access rights for the non-public tables
    in the DES Oracle database.
    
    """
    
    return DESRetriever(RADeg, decDeg, DR = 'Y3', halfBoxSizeDeg = halfBoxSizeDeg, optionsDict = optionsDict)

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

    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    try:
        connection=optionsDict['connection']
    except:
        raise Exception, "No DES database connection"
    
    if decDeg > 5:
        print "... outside DES area - skipping ..."
        return None
    
    outFileName=cacheDir+os.path.sep+"DES%s_%.4f_%.4f_%.4f.csv" % (DR, RADeg, decDeg, halfBoxSizeDeg)      

    # This bit has to be outside, in case we load data from cache
    if DR == 'Y3':
        magKey="SOF_CM_MAG_CORRECTED_$BAND"
        magErrKey="SOF_CM_MAG_ERR_$BAND"
    elif DR == 'DR1':
        magKey="MAG_AUTO_$BAND_DERED"
        magErrKey="MAGERR_AUTO_$BAND"
    else:
        raise Exception, "didn't recognise requested DES DR"

    if os.path.exists(outFileName) == False or 'refetch' in optionsDict.keys() and optionsDict['refetch'] == True:
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        # Y3 Gold v1.0
        #query="SELECT coadd_object_id, ra, dec, ngmix_cm_mag_g - (3.186 * EBV_SFD98) AS cm_mag_g, ngmix_cm_mag_r - (2.140 * EBV_SFD98) AS cm_mag_r, ngmix_cm_mag_i - (1.569 * EBV_SFD98) AS cm_mag_i, ngmix_cm_mag_z - (1.196 * EBV_SFD98) AS cm_mag_z, ngmix_cm_mag_err_g AS cm_mag_err_g, ngmix_cm_mag_err_r AS cm_mag_err_r, ngmix_cm_mag_err_i AS cm_mag_err_i, ngmix_cm_mag_err_z AS cm_mag_err_z FROM y3_gold_1_0 SAMPLE WHERE ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f AND flag_gold = 0 AND flag_footprint = 1 AND flag_foreground = 0 AND extended_class_mash BETWEEN 3 AND 4" % (RAMin, RAMax, decMin, decMax) 
        # Y3 Gold v2.0
        if DR == 'Y3':
            query="SELECT COADD_OBJECT_ID, RA, DEC, DNF_ZMC_SOF, BPZ_ZMC_SOF, SOF_CM_MAG_CORRECTED_G, SOF_CM_MAG_CORRECTED_R, SOF_CM_MAG_CORRECTED_I, SOF_CM_MAG_CORRECTED_Z, SOF_CM_MAG_ERR_G, SOF_CM_MAG_ERR_R, SOF_CM_MAG_ERR_I, SOF_CM_MAG_ERR_Z FROM Y3_GOLD_2_0 WHERE FLAGS_FOOTPRINT = 1 and FLAGS_FOREGROUND = 0 and bitand(FLAGS_GOLD, 62) = 0 and EXTENDED_CLASS_MASH_SOF = 3 and SOF_CM_MAG_I between 16 and 24 AND RA BETWEEN %.6f AND %.6f AND DEC BETWEEN %.6f and %.6f" % (RAMin, RAMax, decMin, decMax)
        elif DR == 'DR1':
            query="SELECT RA, DEC, MAG_AUTO_G_DERED, MAG_AUTO_R_DERED, MAG_AUTO_I_DERED, MAG_AUTO_Z_DERED, MAGERR_AUTO_G, MAGERR_AUTO_R, MAGERR_AUTO_I, MAGERR_AUTO_Z FROM DR1_MAIN WHERE WAVG_SPREAD_MODEL_I + 3.0*WAVG_SPREADERR_MODEL_I > 0.005 and WAVG_SPREAD_MODEL_I + 1.0*WAVG_SPREADERR_MODEL_I > 0.003 and WAVG_SPREAD_MODEL_I - 1.0*WAVG_SPREADERR_MODEL_I > 0.001 and WAVG_SPREAD_MODEL_I > -1 and IMAFLAGS_ISO_I = 0 and MAG_AUTO_I < 24 and RA BETWEEN %.6f AND %.6f AND DEC BETWEEN %.6f and %.6f" % (RAMin, RAMax, decMin, decMax)
        else:
            raise Exception, "didn't recognise requested DES DR"
        
        # This should restart the connection if it drops
        if connection.ping() == False:
            print "... DES database connection lost: reconnecting ..."
            import easyaccess as ea
            connection=ea.connect()
        if connection.ping() == True:
            connection.query_and_save(query, outFileName)
        else:
            raise Exception, "DES database connection lost"

    # No output is written if the query returns nothing... so we'll write a blank file to save repeating the query next time
    catalog=[]
    if os.path.exists(outFileName) == False:
        outFile=file(outFileName, "w")
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
def KIDSDR3Retriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves KIDS DR3 photometry, from the ESO Catalogue Facility. This requires an ESO Portal login,
    which should be stored under $HOME/.zCluster/ESOLogin. That file should just contain the username on the
    first line, and the password on the second. If security is a concern, this could be stored in an encrypted
    directory (e.g., KDE's Vaults facility). See the bin/zCluster script for how the login is handled.
    
    """

    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    try:
        br=optionsDict['connection']
    except:
        raise Exception, "No ESO database connection"
    
    # This defines the rough final KIDS area (split into two fields) - DR3 is a fraction of this
    inKIDSRegion=False
    if RADeg > 120 and RADeg < 240 and decDeg > -5 and decDeg < 5:
        inKIDSRegion=True
    if RADeg > 330 and decDeg > -37 and decDeg < -25:
        inKIDSRegion=True
    if RADeg < 60 and decDeg > -37 and decDeg < -25:
        inKIDSRegion=True
    if inKIDSRegion == False:
        print "... outside KIDS area - skipping ..."
        return None

    outFileName=cacheDir+os.path.sep+"KIDSDR3_%.4f_%.4f_%.4f.fits" % (RADeg, decDeg, halfBoxSizeDeg)      
    if os.path.exists(outFileName) == False or 'refetch' in optionsDict.keys() and optionsDict['refetch'] == True:
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        print "... downloading catalog %s ..." % (outFileName)
        queryURL="https://www.eso.org/qi/catalogQuery/index/112"
        br.open(queryURL)
        br.select_form(nr=1)
        br.form['target']='%.6f %.6f' % (RADeg, decDeg)
        br.form['radius']='%.3f' % (halfBoxSizeDeg)     # default unit is degrees
        br.form['param_112_7']='=0'                     # SG2DPHOT: 0 = galaxy
        br.form.find_control('selectedFields_112').readonly=False
        br.form['selectedFields_112']=',row_112_1,row_112_2,row_112_5,row_112_6,row_112_7,row_112_126,row_112_127,row_112_128,row_112_129,row_112_134,row_112_135,row_112_136,row_112_137,row_112_161,row_112_162,row_112_163,row_112_164,row_112_165,row_112_166,row_112_167,row_112_168,row_112_173,row_112_174,row_112_175,row_112_176,row_112_177,row_112_178'    
        a=br.submit(id = 'exportGeoButtonSpatial')
        b=a.get_data()
        outFile=file(outFileName, "wb")
        outFile.write(b)
        outFile.close()

    # Load/parse table
    print "... reading catalog %s ..." % (outFileName)
    try:
        tab=atpy.Table().read(outFileName)
    except ValueError:
        print "... no objects in catalog - skipping ..."
        return None
    
    magKey="MAG_GAAP_$BAND"
    magErrKey="MAGERR_GAAP_$BAND"
    
    # First, get rid of nans or nonsensical values
    bands=['U', 'G', 'R', 'I']
    for b in bands:
        #tab=tab[np.where(np.isnan(tab[magKey.replace("$BAND", b)]) == False)]
        #tab=tab[np.where(np.isnan(tab[magErrKey.replace("$BAND", b)]) == False)]
        tab=tab[np.where(tab[magKey.replace("$BAND", b)] > 0)]
        tab=tab[np.where(tab[magErrKey.replace("$BAND", b)] > 0)]
        
    idCount=0
    catalog=[]
    for row in tab:
        idCount=idCount+1
        photDict={}
        photDict['id']=idCount    # just so we have something - we could use ID but skipping for now
        photDict['RADeg']=row['RAJ2000']
        photDict['decDeg']=row['DECJ2000']
        photDict['u']=row[magKey.replace("$BAND", "U")]-row['EXT_SFD_U']
        photDict['g']=row[magKey.replace("$BAND", "G")]-row['EXT_SFD_G']
        photDict['r']=row[magKey.replace("$BAND", "R")]-row['EXT_SFD_R']
        photDict['i']=row[magKey.replace("$BAND", "I")]-row['EXT_SFD_I']
        photDict['uErr']=row[magErrKey.replace("$BAND", "U")]
        photDict['gErr']=row[magErrKey.replace("$BAND", "G")]
        photDict['rErr']=row[magErrKey.replace("$BAND", "R")]
        photDict['iErr']=row[magErrKey.replace("$BAND", "I")]
        # Apply mag error cuts if given
        # For PS1, missing values are -999 - our current checkMagErrors routine will fish those out
        # We're just making the mag unconstrained here (missing data), rather than applying a limit
        # If we don't have a minimum of three useful bands, reject
        if 'maxMagError' in optionsDict.keys():
            keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['u', 'g', 'r', 'i'])
        else:
            keep=True
                    
        if keep == True:
            catalog.append(photDict)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def ATLASDR3Retriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves VST ATLAS DR3 photometry, from the ESO Catalogue Facility. This requires an ESO Portal login,
    which should be stored under $HOME/.zCluster/ESOLogin. That file should just contain the username on the
    first line, and the password on the second. If security is a concern, this could be stored in an encrypted
    directory (e.g., KDE's Vaults facility). See the bin/zCluster script for how the login is handled.
    
    """

    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    try:
        br=optionsDict['connection']
    except:
        raise Exception, "No ESO database connection"
    
    # This defines the rough ATLAS DR3 area
    inATLASRegion=False
    if RADeg > 21*15 or RADeg < 4.2*60:
        if decDeg > -42 and decDeg < -7:
            inATLASRegion=True
    if inATLASRegion == False:
        print "... outside ATLAS area - skipping ..."
        return None

    outFileName=cacheDir+os.path.sep+"ATLASDR3_%.4f_%.4f_%.4f.fits" % (RADeg, decDeg, halfBoxSizeDeg)      
    if os.path.exists(outFileName) == False or 'refetch' in optionsDict.keys() and optionsDict['refetch'] == True:
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        print "... downloading catalog %s ..." % (outFileName)
        queryURL="https://www.eso.org/qi/catalogQuery/index/147"
        br.open(queryURL)
        br.select_form(nr=1)
        br.form['target']='%.6f %.6f' % (RADeg, decDeg)
        br.form['radius']='%.3f' % (halfBoxSizeDeg)     # default unit is degrees
        #br.form['param_112_7']='=0'                     # SG2DPHOT: 0 = galaxy
        br.form.find_control('selectedFields_147').readonly=False
        # Selected Petrosian mags below (and extinction corrections needed)
        br.form['selectedFields_147']=',row_147_5,row_147_6,row_147_30,row_147_31,row_147_35,row_147_36,row_147_37,row_147_38,row_147_39,row_147_41,row_147_42,row_147_65,row_147_66,row_147_89,row_147_90,row_147_113,row_147_114,row_147_137,row_147_138'    
        a=br.submit(id = 'exportGeoButtonSpatial')
        b=a.get_data()
        outFile=file(outFileName, "wb")
        outFile.write(b)
        outFile.close()

    # Load/parse table
    print "... reading catalog %s ..." % (outFileName)
    try:
        tab=atpy.Table().read(outFileName)
    except ValueError:
        print "... no objects in catalog - skipping ..."
        return None
    
    magKey="$BANDPETROMAG"
    magErrKey="$BANDPETROMAGERR"
    
    # First, get rid of nans or nonsensical values
    bands=['U', 'G', 'R', 'I', 'Z']
    for b in bands:
        tab[magKey.replace("$BAND", b)][np.isnan(tab[magKey.replace("$BAND", b)])]=99.
        tab[magErrKey.replace("$BAND", b)][np.isnan(tab[magErrKey.replace("$BAND", b)])]=99.

    # Star-galaxy separation (could be done at query level...)
    tab=tab[np.where(tab['PGALAXY'] > 0.5)]

    idCount=0
    catalog=[]
    for row in tab:
        idCount=idCount+1
        photDict={}
        photDict['id']=idCount    # just so we have something - we could use ID but skipping for now
        photDict['RADeg']=row['RA2000']
        photDict['decDeg']=row['DEC2000']
        photDict['u']=row[magKey.replace("$BAND", "U")]-row['AU']
        photDict['g']=row[magKey.replace("$BAND", "G")]-row['AG']
        photDict['r']=row[magKey.replace("$BAND", "R")]-row['AR']
        photDict['i']=row[magKey.replace("$BAND", "I")]-row['AI']
        photDict['z']=row[magKey.replace("$BAND", "Z")]-row['AZ']
        photDict['uErr']=row[magErrKey.replace("$BAND", "U")]
        photDict['gErr']=row[magErrKey.replace("$BAND", "G")]
        photDict['rErr']=row[magErrKey.replace("$BAND", "R")]
        photDict['iErr']=row[magErrKey.replace("$BAND", "I")]
        photDict['zErr']=row[magErrKey.replace("$BAND", "Z")]
        # Apply mag error cuts if given
        # For PS1, missing values are -999 - our current checkMagErrors routine will fish those out
        # We're just making the mag unconstrained here (missing data), rather than applying a limit
        # If we don't have a minimum of three useful bands, reject
        if 'maxMagError' in optionsDict.keys():
            keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['u', 'g', 'r', 'i', 'z'])
        else:
            keep=True
                    
        if keep == True:
            catalog.append(photDict)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def PS1Retriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves PS1 photometry at the given position.
        
    """

    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)

    outFileName=cacheDir+os.path.sep+"PS1_%.4f_%.4f_%.4f.xml" % (RADeg, decDeg, halfBoxSizeDeg)      
    print "... getting PS1 photometry (file: %s) ..." % (outFileName)

    if decDeg < -30:
        print "... outside PS1 area - skipping ..."
        return None
    
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
        #tab=tab[np.where(tab['rMeanPSFMag']-tab['rMeanKronMag'] >= 0.5)]
        for row in tab:
            idCount=idCount+1
            photDict={}
            photDict['id']=idCount    # just so we have something - we could use PS1 ID but skipping for now
            photDict['RADeg']=row['raMean']
            photDict['decDeg']=row['decMean']
            photDict['g']=row['gMeanKronMag']
            photDict['r']=row['rMeanKronMag']
            photDict['i']=row['iMeanKronMag']
            photDict['z']=row['zMeanKronMag']
            photDict['Y']=row['yMeanKronMag']
            photDict['gErr']=row['gMeanKronMagErr']
            photDict['rErr']=row['rMeanKronMagErr']
            photDict['iErr']=row['iMeanKronMagErr']
            photDict['zErr']=row['zMeanKronMagErr']
            photDict['YErr']=row['yMeanKronMagErr']
            #photDict['g']=row['gMeanApMag']
            #photDict['r']=row['rMeanApMag']
            #photDict['i']=row['iMeanApMag']
            #photDict['z']=row['zMeanApMag']
            #photDict['Y']=row['yMeanApMag']
            #photDict['gErr']=row['gMeanApMagErr']
            #photDict['rErr']=row['rMeanApMagErr']
            #photDict['iErr']=row['iMeanApMagErr']
            #photDict['zErr']=row['zMeanApMagErr']
            #photDict['YErr']=row['yMeanApMagErr']
            
            # Correct for dust extinction
            # Taken from: http://www.mso.anu.edu.au/~brad/filters.html
            photDict['g']=photDict['g']-EBMinusV*3.322
            photDict['r']=photDict['r']-EBMinusV*2.544 
            photDict['i']=photDict['i']-EBMinusV*2.265 
            photDict['z']=photDict['z']-EBMinusV*1.846 
            photDict['Y']=photDict['Y']-EBMinusV*1.570 

            # Apply mag error cuts if given
            # For PS1, missing values are -999 - our current checkMagErrors routine will fish those out
            # We're just making the mag unconstrained here (missing data), rather than applying a limit
            # If we don't have a minimum of three useful bands, reject
            if 'maxMagError' in optionsDict.keys():
                keep=checkMagErrors(photDict, optionsDict['maxMagError'], bands = ['g', 'r', 'i', 'z'])
            else:
                keep=True
            
            # Additional PS1 quality cut
            # NOTE: this does not help (makes scatter worse actually)
            #if keep == True and row['gQfPerfect'] > 0.9 and row['rQfPerfect'] > 0.9 and row['iQfPerfect'] > 0.9 and row['gQfPerfect'] > 0.9:
                #keep=True
            #else:
                #keep=False
            
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
        
    if os.path.exists(outFileName) == False or 'refetch' in optionsDict.keys() and optionsDict['refetch'] == True:
        
        print "... fetching from the internet ..."
        
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
                    print "... problem with file %s - removing and retrying in 70 sec ..." % (outFileName)
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
                if 'maxMagError' in optionsDict.keys():
                    keep=checkMagErrors(photDict, optionsDict['maxMagError'])
                else:
                    keep=True
                if keep == True:
                    catalog.append(photDict)
        
    return catalog

#-------------------------------------------------------------------------------------------------------------
def DECaLSRetriever(RADeg, decDeg, halfBoxSizeDeg = 18.0/60.0, optionsDict = {}):
    """Retrieves DECaLS DR5 tractor catalogs (if they exist) at the given position
    
    NOTE: we're using some DR6 files here... but DR5 catalogs... it seems like DR6 is missing the actual 
    DECaLS data?
    
    """

    makeCacheDir()
    
    if 'altCacheDir' in optionsDict.keys():
        cacheDir=optionsDict['altCacheDir']
    else:
        cacheDir=CACHE_DIR
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)

    # Organised such that after this, have subdir with degrees RA (e.g. 000/ is 0 < RADeg < 1 degree)
    basePath="http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr5/tractor/"

    outFileName=cacheDir+os.path.sep+"DECaLSDR5_%.4f_%.4f_%.2f.fits" % (RADeg, decDeg,
                                                                            halfBoxSizeDeg)
    
    print "... getting DECaLS photometry (file: %s) ..." % (outFileName)

    bricksTab=optionsDict['bricksTab']
    DR6Tab=optionsDict['DR6Tab']
    
    if os.path.exists(outFileName) == False:
        
        # Find matching tractor catalogs and download (temporary) .fits table files
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        mask=np.logical_and(np.greater(bricksTab['RA1'], RAMin), np.less(bricksTab['RA2'], RAMax))
        mask=np.logical_and(mask, np.greater(bricksTab['DEC1'], decMin))
        mask=np.logical_and(mask, np.less(bricksTab['DEC2'], decMax))
        matchTab=bricksTab[np.where(mask)]
        count=0
        tractorTabs=[]
        cacheFileNames=[]
        matchTab=atpy.join(matchTab, DR6Tab, keys = 'BRICKNAME')
        for mrow in matchTab:
            print "... retrieving tractor catalog from web ..."
            url=basePath+"%03d" % np.floor(mrow['RA'])
            url=url+os.path.sep+"tractor-"+mrow['BRICKNAME']+".fits"
            fileName=cacheDir+os.path.sep+"tractor-tmp_%s.fits" % (mrow['BRICKNAME'])
            cacheFileNames.append(fileName)
            urllib.urlretrieve(url, filename = fileName)
            try:
                tractorTab=atpy.Table.read(fileName)
                tractorTabs.append(tractorTab)
            except:
                print "... probably a 404 error for %s ..." % (fileName)
        # Stitch catalogs together and write to cache, delete temporary files
        if len(tractorTabs) > 0:
            uberTab=atpy.vstack(tractorTabs)
            uberTab.write(outFileName)
            for fileName in cacheFileNames:
                os.remove(fileName)
    
    # Load catalog
    if os.path.exists(outFileName) == True:
        
        tab=atpy.Table().read(outFileName)
        
        # Remove stars
        tab=tab[np.where(tab['type'] != 'PSF')]
        
        # WISE fluxes are available...
        bands=['u', 'g', 'r', 'i', 'z', "w1", "w2"]# , 'Y']

        # Convert nanomaggies to mags and do extinction correction
        bricksInTab=np.unique(tab['brickname'])
        for brickName in bricksInTab:
            brickExtTab=DR6Tab[np.where(DR6Tab['BRICKNAME'] == brickName)]
            brickIndices=np.where(tab['brickname'] == brickName)
            brickMask=np.zeros(len(tab), dtype = bool)
            brickMask[brickIndices]=True
            for b in bands:
                magLabel='%s_mag' % (b)
                magErrLabel='%s_magErr' % (b)
                extKey='ext_%s' % (b)
                colsToAdd=[magLabel, magErrLabel]
                for col in colsToAdd:
                    if col not in tab.keys():
                        tab.add_column(atpy.Column(np.ones(len(tab))*99., col))
                nanoMaggies=np.array(tab['flux_%s' % (b)])
                validMask=np.greater(nanoMaggies, 0.)
                mask=np.logical_and(brickMask, validMask)
                tab[magLabel][mask]=-2.5*np.log10(nanoMaggies[mask])+22.5
                validFluxErr=np.sqrt(1./tab['flux_ivar_%s' % (b)][mask])
                tab[magErrLabel][mask]=1./(nanoMaggies[mask]/validFluxErr)
                if extKey in brickExtTab.keys():
                    dustCorrMag=brickExtTab[extKey][0]
                    tab[magLabel][mask]=tab[magLabel][mask]-dustCorrMag
        
        catalog=[]
        for row in tab: 
            photDict={}
            photDict['id']=row['objid']
            photDict['RADeg']=row['ra']
            photDict['decDeg']=row['dec']
            for b in bands:
                photDict[b]=row['%s_mag' % (b)]
                photDict[b+"Err"]=row['%s_magErr' % (b)]
            if 'maxMagError' in optionsDict.keys():
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
                    print "... problem with file %s - removing and retrying in 70 sec ..." % (outFileName)
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
