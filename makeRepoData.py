import numpy as np
#from sunpy.time import parse_time
import pickle
from makeMegaDataStructure import AllRes


f =  open('allRes.pkl', 'rb')
res, id2time, time2time = pickle.load(f)
f.close()  
keys = np.sort(np.array([key for key in res.keys()]))

# Projected CORSET values
if False:
    f1 = open('wombatProjValues.dat', 'w')
    f1.write('           Time        CORIDA        CORIDB      MA      hA    CPAA      AWA     velA      MB      hB    CPAB      AWB     velB \n')
    for key in keys:
        aRes = res[key]
    
        if aRes.CORSETidA or aRes.CORSETidB:
            idA = str(aRes.CORSETidA).rjust(14)
            idB = str(aRes.CORSETidB).rjust(14)
        
            outstr = key + idA + idB
        
        
            if aRes.CORSETmassMA:
                massA = str(aRes.CORSETmassMA)
                idx = massA.find('.')
                if (len(massA) - idx) == 2:
                    massA = massA+'0'
                massA = massA.rjust(8)               
                htA = str(aRes.CORSEThtMA).rjust(8)
                cpaA = str(int(aRes.CORSETcpaMA)).rjust(8)
                angA = str(int(aRes.CORSETangMA)).rjust(9)
                velA = str(int(aRes.CORSETvelMA)).rjust(9)
                outstr = outstr + massA + htA + cpaA + angA + velA
            else:
                outstr = outstr +'    None    None    None     None     None'
            
   
            if aRes.CORSETmassMB:
                massB = str(aRes.CORSETmassMB)
                idx = massB.find('.')
                if (len(massB) - idx) == 2:
                    massB = massB+'0'
                massB = massB.rjust(8)
                htB = str(aRes.CORSEThtMB).rjust(8)
                cpaB = str(int(aRes.CORSETcpaMB)).rjust(8)
                angB = str(int(aRes.CORSETangMB)).rjust(9)
                velB = str(int(aRes.CORSETvelMB)).rjust(9)
                outstr = outstr + massB + htB + cpaB + angB + velB
            else:
                outstr = outstr +'    None    None    None     None     None'
            f1.write(outstr+'\n')
    f1.close()

# Deprojected CORSET values
if False:
    f1 = open('wombatDeprojValues.dat', 'w')
    f1.write('           Time        CORIDA        CORIDB   DP1_MA   DP1_MB DP1_Lat DP1_Lon        DP2_M DP2_LatA DP2_LatB DP2_Lon        DP3_MA   DP3_MB DP3_LatA DP3_LatB  DP3_Lon\n')
    for key in keys:
        aRes = res[key]
        if aRes.deprojTimes[0]:
            idA = str(aRes.CORSETidA).rjust(14)
            idB = str(aRes.CORSETidB).rjust(14)
            outstr = key + idA + idB
            
            # DP1 (LB)
            if aRes.LBdeprojLat:
                outstr = outstr +'{:9.2f}'.format(aRes.LBdeprojMA) +'{:9.1f}'.format(aRes.LBdeprojMB) + '{:8.2f}'.format(aRes.LBdeprojLat) + '{:8.1f}'.format(aRes.LBdeprojLon)
            else: 
                 outstr = outstr +'     None     None    None    None'
                 
            
            # DP2
            maxM = np.max(aRes.deprojMasses)
            maxIdx = np.where(aRes.deprojMasses == maxM)[0]
            maxM = '{:13.2f}'.format(maxM)
            latA = '{:9.1f}'.format(float(aRes.deprojLatA[maxIdx[0]]))
            latB = '{:9.1f}'.format(float(aRes.deprojLatB[maxIdx[0]]))
            lon = '{:9.1f}'.format(float(aRes.deprojLons[maxIdx[0]]))
            outstr = outstr + maxM + latA + latB + lon
            
            # DP3
            maxMA = np.max(aRes.CdeprojMassesA)
            maxMB = np.max(aRes.CdeprojMassesB)
            maxIdxA = np.where(aRes.CdeprojMassesA == maxMA)[0]
            maxMA = '{:13.2f}'.format(maxMA)
            maxIdxB = np.where(aRes.CdeprojMassesB == maxMB)[0]
            maxMB = '{:9.2f}'.format(maxMB)
            latA = '{:9.1f}'.format(float(aRes.CdeprojLatA[maxIdxA[0]]))
            latB = '{:9.1f}'.format(float(aRes.CdeprojLatB[maxIdxB[0]]))
            lon = '{:9.1f}'.format(float(aRes.CdeprojLons[maxIdx[0]]))
            outstr = outstr +  maxMA + maxMB + latA + latB + lon
            
            f1.write(outstr+'\n')
    f1.close()
    
# GCS values
if False:
    f1 = open('wombatGCSValues.dat', 'w')
    f1.write('           Time        CORIDA        CORIDB   SourceCat        SourceTime     Lat     Lon    Tilt      AW   kappa       h     vel       MA       MB\n')
    for key in keys:
        aRes = res[key]
        allNames = []
        if aRes.GCSnamesA:
            allNames.append(aRes.GCSnamesA) 
        if aRes.GCSnamesB:
            allNames.append(aRes.GCSnamesB) 
        allNames = np.unique(allNames)
        for name in allNames:
            idA = str(aRes.CORSETidA).rjust(14)
            idB = str(aRes.CORSETidB).rjust(14)
            outstr = key + idA + idB + name.rjust(12)
            
            idxA, idxB = None, None
            GCStime = None
            if aRes.GCSnamesA:
                if name in aRes.GCSnamesA:
                    idxA = np.where(np.array(aRes.GCSnamesA) == name)[0]
                    GCStime = aRes.GCStimesA[idxA[0]]
            if aRes.GCSnamesB:
                if name in aRes.GCSnamesB:
                    idxB = np.where(np.array(aRes.GCSnamesB) == name)[0]
                    GCStime = aRes.GCStimesB[idxB[0]]
            
            outstr = outstr + '   ' + GCStime
            
            GCSvals = aRes.GCSvals[name]
            GCSstr = GCSvals[0].rjust(8) +  GCSvals[1].rjust(8) +  GCSvals[2].rjust(8) +  GCSvals[3].rjust(8) +  GCSvals[4].rjust(8) +  GCSvals[5].rjust(8) +  GCSvals[6].rjust(8)
            outstr = outstr + GCSstr
            
            if type(idxA) != type(None):
                outstr = outstr + '{:9.2f}'.format(aRes.GCSmassesA[idxA[0]])
            else:
                outstr = outstr +'     None'

            if type(idxB) != type(None):
                outstr = outstr + '{:9.2f}'.format(aRes.GCSmassesB[idxB[0]])
            else:
                outstr = outstr +'     None'
            
            #print (outstr)
            f1.write(outstr+'\n')
    f1.close()
    