import numpy as np
#from sunpy.time import parse_time
import pickle
from makeMegaDataStructure import AllRes
import os


f =  open('allRes.pkl', 'rb')
res, id2time, time2time = pickle.load(f)
f.close()  
keys = np.sort(np.array([key for key in res.keys()]))

pngPath = '/Users/kaycd1/STEREO_Mass/massProfFigs/'
txtPath = '/Users/kaycd1/STEREO_Mass/indivFiles/'

yrs = ['2007','2008','2009','2010','2011','2012','2013','2014']
for yr in yrs:
    f = np.genfromtxt('corset_cme_list_'+yr+'.csv', delimiter=',', dtype='str')
    f2 = open('madcat_cme_list_'+yr+'.csv', 'w')
    nPts = len(f[:,0])
    nThings = len(f[0,:])
    for i in range(nPts):
        corID = f[i,1].replace('c','')
        pMass = 'None'
        dMass = 'None'
        lat   =  'None'
        lon   =  'None'
        myPng = 'None'
        myTxt = 'None'
        if corID in id2time.keys():
            key = id2time[corID]
            myRes = res[key]
            if 'A' in corID:
                if myRes.CORSETmassMA:
                    pMass = '{:.2f}'.format(myRes.CORSETmassMA)
                if myRes.CdeprojMassesA[0]:
                    dMass = np.max(myRes.CdeprojMassesA)
                    idxM  = np.where(np.array(myRes.CdeprojMassesA) == dMass)[0]
                    if dMass:
                        dMass = '{:.2f}'.format(dMass)
                    lat   = '{:.2f}'.format(myRes.CdeprojLatA[idxM[0]])
                    lon   = '{:.2f}'.format(myRes.CdeprojLons[idxM[0]])
            if 'B' in corID:
                if myRes.CORSETmassMB:
                    pMass = '{:.2f}'.format(myRes.CORSETmassMB)
                if myRes.CdeprojMassesB[0]:
                    dMass = np.max(myRes.CdeprojMassesB)
                    idxM  = np.where(np.array(myRes.CdeprojMassesB) == dMass)[0]
                    if dMass:
                        dMass = '{:.2f}'.format(dMass)
                    lat   = '{:.2f}'.format(myRes.CdeprojLatB[idxM[0]])
                    lon   = '{:.2f}'.format(myRes.CdeprojLons[idxM[0]])
                
            # add in the file names
            if os.path.exists(pngPath+'massProfile_'+key+'.png'): 
                myPng = 'massProfile_'+key+'.png'
            if os.path.exists(txtPath+key+'.txt'): 
                myTxt = key+'.txt'
        
        outs = ''
    
        for j in range(nThings):
            outs = outs + f[i,j] + ','   
        if i == 0:     
            outs = outs + 'proj Mass,deproj Mass,lat,lon,pngFile,txtFile'
        else:    
            outs = outs + pMass + ',' + dMass  + ',' +  lat  + ',' +  lon + ',' +  myPng  + ',' +  myTxt
        #print(outs)
        f2.write(outs+'\n')
 