""" 

This module contains tools for handling catalogs.
    
"""

from astLib import *
import astropy.table as atpy
import numpy as np
import operator
import os
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import sys
import time

# For adding meta data to output
import datetime

#------------------------------------------------------------------------------------------------------------
def catalog2DS9(catalog, outFileName, constraintsList = [], addInfo = [], idKeyToUse = 'name', 
                RAKeyToUse = 'RADeg', decKeyToUse = 'decDeg', color = "cyan"):
    """Converts a catalog containing object dictionaries into a ds9 region file. Objects will be labelled
    in the .reg file according to the idKeyToUse.
    
    If color == 'key', then use 'color' key in object dictionary to set color.
    
    constraintsList works the same way as selectFromCatalog function
    
    """

    # Cut catalog according to constraints
    cutCatalog=selectFromCatalog(catalog, constraintsList) 
    
    outFile=open(outFileName, "w")
    timeStamp=datetime.datetime.today().date().isoformat()
    comment="# DS9 region file\n"
    outFile.write(comment)
    outFile.write('global dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    for obj in cutCatalog:
        if len(addInfo) > 0:
            infoString=""
            for d in addInfo:
                if infoString != "":
                    infoString=infoString+" "
                if obj[d['key']] != None:
                    infoString=infoString+d['fmt'] % (obj[d['key']])
                else:
                    infoString=infoString+"%s" % (str(obj[d['key']]))
            infoString=" ["+infoString+"]"
        else:
            infoString=""
        if color == 'key':
            colorString=obj['color']
        else:
            colorString=color
        outFile.write("fk5;point(%.6f,%.6f) # point=boxcircle color={%s} text={%s%s}\n" \
                    % (obj[RAKeyToUse], obj[decKeyToUse], colorString, obj[idKeyToUse], infoString))
    outFile.close()
        
#------------------------------------------------------------------------------------------------------------
def selectFromCatalog(catalog, constraintsList):
    """Given a catalog (list of dictionaries representing objects), return a list of objects matching the
    given constraintsList. Each item in constraintsList is a string in the form:
    
    "key < value", "key > value", etc.
    
    Note that the spaces between key, operator (e.g. '<') and value are essential!
    
    """
    
    passedConstraint=catalog
    for constraintString in constraintsList:
        lastCatalog=passedConstraint
        passedConstraint=[]
        for obj in lastCatalog:         
            key, op, value=constraintString.split()
            if eval("obj['%s'] %s %s" % (key, op, value)) == True:
                passedConstraint.append(obj)
    
    return passedConstraint

#------------------------------------------------------------------------------------------------------------
def writeCatalog(catalog, outFileName, keysToWrite, keyFormats, constraintsList, headings = True, 
                 extraHeaderText = None):
    """Dumps the catalog to a .csv.
    
    constraintsList works as in the selectFromCatalog function.
        
    """
    
    # Cut catalog according to constraints
    cutCatalog=selectFromCatalog(catalog, constraintsList)                                           
    
    outFile=open(outFileName, "w")
    
    # Add meta data
    timeStamp=datetime.datetime.today().date().isoformat()
    outFile.write("# Date generated: %s\n" % (timeStamp))
    if extraHeaderText != None:
        outFile.write(extraHeaderText)
    if headings == True:
        heading="# "
        for k in keysToWrite:
            if heading != "# ":
                heading=heading+"\t"
            heading=heading+k
        outFile.write(heading+"\tnotes\n")
        
    count=0
    for obj in cutCatalog:
        count=count+1
        obj['id']=count
        line=""
        for k, f in zip(keysToWrite, keyFormats):
            if line != "":
                line=line+"\t"
            if type(obj[k]) == list:
                if obj[k][0] != None:   # merged cat, just take first item for now
                    line=line+f % (obj[k][0])   
                else:
                    line=line+str(None)
            else:
                if obj[k] != None:
                    try:
                        line=line+f % (obj[k]) 
                    except:
                        print("Argh!")
                        ipshell()
                        sys.exit()
                else:
                    line=line+str(None)
        # Add on a 'notes' column - any key which is just a bool gets added to a , delimited list if True
        notes=""
        for key in list(obj.keys()):
            if type(obj[key]) == bool and obj[key] == True:
                if notes != "":
                    notes=notes+","
                notes=notes+key
        outFile.write(line+"\t"+notes+"\n")
    outFile.close()    

#-------------------------------------------------------------------------------------------------------------
def writeRedshiftsCatalog(catalog, outFileName):
    """Writes a .fits table with cluster photo-zs in it
    
    """
    
    tab=atpy.Table()
    keys=['name', 'RADeg', 'decDeg', 'origRADeg', 'origDecDeg', 'offsetArcmin', 'offsetMpc', 'z', 'delta', 'errDelta', 'CS', 'A']
    for key in keys:
        arr=[]
        for obj in catalog:
            arr.append(obj[key])
        tab.add_column(atpy.Column(arr, key))
    
    # Since generally we cross match against sourcery tables, we may as well zap the non-results as we don't need
    # The try... except bit here is because we've also used this routine to write the nullTest catalog, which doesn't have zs
    try:
        tab=tab[np.where(tab['z'] > 0.)]
    except:
        pass
    
    tab.table_name="zCluster"
    tab.write(outFileName, overwrite = True)
    
#-------------------------------------------------------------------------------------------------------------
def getNEDInfo(obj, nedDir = "NEDResults", radiusDeg = 5.0/60.0, crossMatchRadiusDeg = 2.5/60, refetch = True):
    """Queries NED for matches near each obj (must have keys name, RADeg, decDeg) and returns the nearest 
    cluster match and its distance in arcmin. We search a box of radiusDeg around the object.
        
    """
        
    if os.path.exists(nedDir) == False:
        os.makedirs(nedDir, exist_ok = True)
                    
    name=obj['name']
    RADeg=obj['RADeg']
    decDeg=obj['decDeg']
            
    outFileName=nedDir+os.path.sep+name.replace(" ", "_")+"_%.1f.txt" % (radiusDeg*60.0)        
    if os.path.exists(outFileName) == False or refetch == True:
        print("... fetching NED info for %s ..." % (name))
        try:                
            urllib.request.urlretrieve("http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.6fd&lat=%.6fd&radius=%.2f&dot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&in_objtypes1=QSO&in_objtypes2=Radio&in_objtypes2=SmmS&in_objtypes2=Infrared&in_objtypes2=Xray&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=ascii_tab&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (RADeg, decDeg, radiusDeg*60.0), filename = outFileName)
        except:
            print("WARNING: couldn't get NED info")
            #IPython.embed()
            #sys.exit()
            outFileName=None
            
    nedObjs=parseNEDResult(outFileName)
    
    # Flag matches against clusters - choose nearest one
    rMin=10000
    clusterMatch={}
    if len(nedObjs['RAs']) > 0:
        obj['NED_allClusterMatches']=[]
        for i in range(len(nedObjs['RAs'])):
            ned=nedObjs
            if ned['sourceTypes'][i] == 'GClstr':
                r=astCoords.calcAngSepDeg(ned['RAs'][i], ned['decs'][i], RADeg, decDeg)
                if r < crossMatchRadiusDeg:
                    obj['NED_allClusterMatches'].append(ned['names'][i])
                if r < rMin and r < crossMatchRadiusDeg:
                    keepName=False
                    if 'name' in clusterMatch:
                        if "ABELL" in clusterMatch['name']:
                            keepName=True
                    if keepName == False:
                        rMin=r
                        clusterMatch['name']=ned['names'][i]
                        if ned['redshifts'][i] != 'N/A':
                            clusterMatch['z']=float(ned['redshifts'][i])
                        else:
                            clusterMatch['z']=None
                        clusterMatch['rArcmin']=rMin*60.0
                        clusterMatch['NED_RADeg']=float(ned['RAs'][i])
                        clusterMatch['NED_decDeg']=float(ned['decs'][i])
    
    if clusterMatch != {}:
        return clusterMatch
    else:
        return None

#-------------------------------------------------------------------------------------------------------------
def parseNEDResult(inFileName, onlyObjTypes = None):
    """Parses NED tab-delimited text file query result, returns dictionary.
    
    onlyObjTypes can be a string indicating types of objects only to include e.g. GClstr
    
    """
    
    if inFileName != None and os.path.exists(inFileName):
        inFile=open(inFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    else:
        # Fail safe in case we couldn't contact NED
        lines=[]

    dataStarted=False
    labels=[]
    names=[]
    RAs=[]
    decs=[]
    sourceTypes=[]
    redshifts=[]
    for line in lines:
        bits=line.split("\t")
        if bits[0] == "1":
            dataStarted=True
        if dataStarted == True:
            if onlyObjTypes == str(bits[4]) or onlyObjTypes == None:
                labels.append(bits[0])
                names.append(bits[1])
                RAs.append(float(bits[2]))
                decs.append(float(bits[3]))
                sourceTypes.append(str(bits[4]))
                if bits[6] == '':
                    redshifts.append('N/A')
                else:
                    redshifts.append(str(bits[6]))

    return {'labels': labels, 'names': names, 'RAs': RAs, 'decs': decs, 'sourceTypes': sourceTypes, 'redshifts': redshifts}

#-------------------------------------------------------------------------------------------------------------
def makeRADecString(RADeg, decDeg):
    """Switched to using %.5f_%.5f as part of image file names.
    
    """
    return "%.5f_%.5f" % (RADeg, decDeg)
