import numpy as np
from sunpy.coordinates import HeliocentricEarthEcliptic, get_body_heliographic_stonyhurst, get_horizons_coord
from sunpy.time import parse_time
import pickle
from makeMegaDataStructure import AllRes

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)



def elTheory(Rin,theta,limb=0.63,center=False, returnAll=False):
    # units are in solar radii and degrees (based on IDL)
    radeg = 180. / np.pi
    # Port of SSWIDL code
    # 0.63 is default for limb darkening
    u = limb
    
    #const = sigma * pi / 2
    #where sigma is the scattering cross section (7.95e-26 per steradian)
    const = 1.24878e-25
    if not center:
        const = const/(1-u/3) # convert to mean solar brightness
    
    # make sure theta is less than 90 deg
    #if theta >= 90:
    #    print ('PoS angle greater than 90, exiting elTheory without running')
    #    return None,  None,  None,  None,  None
    
    R = Rin/np.cos(theta/radeg)
    sinchi2 = (Rin/R)**2	# angle between sun center, electron, observer
    s = 1./R
    #if s >= 1:
    #    print ('Error in computing s, exiting elTheory without running')
    #    return None,  None,  None,  None,  None
    s2 = s**2
    c2 = (1.-s2)
    #if c2 < 0:
    #    print ('Error in computing c2, exiting elTheory without running')
    #    return None,  None,  None,  None,  None	
    c = np.sqrt(c2)			# cos(omega)
    g = c2*(np.log((1.+s)/c))/s
    
    #  Compute Van de Hulst Coefficients
    #  Expressions are given in Billings (1968) after Minnaert (1930)
    ael = c*s2
    cel = (4.-c*(3.+c2))/3.
    bel = -(1.-3.*s2-g*(1.+3.*s2))/8.
    del0 = (5.+s2-g*(5.-s2))/8.
    
    #  Compute electron brightness
    #  pB is polarized brightness (Bt-Br)
    Bt = const*( cel + u*(del0-cel) )
    pB = const* sinchi2 *( (ael + u*(bel-ael) ) )
    Br = Bt-pB
    B = Bt+Br
    Pol = pB/B
    
    if returnAll:
        return R,B,Bt,Br,Pol
    else:
        return R,B

def findDir(Aparams, Bparams):
    # params are [satLon, h, mass]
    lonA, hA, mA = Aparams
    lonB, hB, mB = Bparams

    # Get plane of sky in 0-360 deg
    posA = (lonA + 90 + 360) % 360
    posB = (lonB + 90 + 360) % 360
    if hA > 0:
        arng = [(posA-90 +360) % 360, (posA+90 +360) % 360]
    else:
        arng = [(posA-270 +360), (posA-90 + 360) % 360]
    if arng[1] < arng[0]: arng[1] += 360
        
    if hB > 0:
        brng = [(posB-90 + 360) % 360, (posB+90 + 360) % 360]
    else:
        brng = [(posB-270 +360) % 360, (posB-90 + 360) % 360]
        
    if brng[1] < brng[0]: brng[1] += 360
    
    lonMin = np.max([arng[0], brng[0]])    
    lonMax = np.min([arng[1], brng[1]])   
    if lonMax < lonMin:
        lonMax += 360 
    
    # Get range to check
    # h + for west/right and - for east/left 
    '''lonMin, lonMax = np.min([lonA, lonB]), np.max([lonA, lonB])
    if (hA > 0) & (hB > 0):
        # higher lon than both
        angMin, angMax = lonMax, lonMin+180
    elif (hA < 0) & (hB < 0):
        # lower lon than both
        angMin, angMax = lonMax-180, lonMin
    elif (hA < 0) & (hB > 0):
        # lower lon A and higher than B
        aPOS = lonA - 90
        bPOS = lonB + 90
        # check if both POS are on backside
        if (np.abs(aPOS) > 90) & (np.abs(bPOS) > 90):
            aPOS = (aPOS+360) % 360 
            bPOS = (bPOS+360) % 360 
        arng = [aPOS-90, aPOS+90]
        brng = [bPOS-90, bPOS+90]
        angMin = np.max([arng[0], brng[0]])
        angMax = np.min([arng[1], brng[1]])
    elif (hA > 0) & (hB < 0):
        # higher lon A and lower than B'
        aPOS = lonA + 90
        bPOS = lonB - 90
        # check if both POS are on backside
        if (np.abs(aPOS) > 90) & (np.abs(bPOS) > 90):
            aPOS = (aPOS+360) % 360 
            bPOS = (bPOS+360) % 360 
        arng = [aPOS-90, aPOS+90]
        brng = [bPOS-90, bPOS+90]
        angMin = np.max([arng[0], brng[0]])
        angMax = np.min([arng[1], brng[1]])'''
    dirs = np.arange(lonMin,lonMax,1.)
    sepA = dirs - (lonA + np.sign(hA)*90)
    sepB = dirs - (lonB + np.sign(hB)*90)
    
    sepA[np.where(sepA > 180.)] = 360 - sepA[np.where(sepA > 180.)] 
    sepB[np.where(sepB > 180.)] = 360 - sepB[np.where(sepB > 180.)] 
    
    routA0, BA0 = elTheory(hA, 0)
    routB0, BB0 = elTheory(hB, 0)
    
    routAs, BAs = elTheory(hA, np.abs(sepA))
    routBs, BBs = elTheory(hB, np.abs(sepB))
    
    normAs = BAs/BA0
    normBs = BBs/BB0

    massA = mA/normAs
    massB = mB/normBs
    massDiff = massA - massB
    

    # get the first cut at finding the direction
    # need to be extra careful on first pass to not pull where
    # it starts blowing up wildly but passes through zero
    minIdx = np.where(np.abs(massDiff) == np.min(np.abs(massDiff)))[0]
    '''deltaMD = np.zeros(len(massDiff))
    deltaMD[1:] = np.abs(massDiff[1:] - massDiff[:-1])
    medDMD  = np.median(deltaMD)
    deltaMD[0] = deltaMD[1]
    niceIdx = np.where(np.abs(deltaMD/medDMD) <= 5)[0]    
    newIdx = np.where(np.abs(massDiff) == np.min(np.abs(massDiff[niceIdx])))[0]'''
    #for i in range(len(sepA)):
    #    print (i, dirs[i], massDiff[i])
        
    newIdx = minIdx
    dir1 = dirs[newIdx[0]]

    # need to add check for at endpoints
    if massDiff[newIdx[0]]*massDiff[newIdx[0]+1] < 1:
        dir2 = dirs[newIdx[0]+1]
    else:
        dir2 = dir1
        dir1 = dirs[newIdx[0]-1]
        
    # Refine the fit
    dirs = np.linspace(dir1,dir2,100)
    sepA = dirs - (lonA + np.sign(hA)*90)
    sepB = dirs - (lonB + np.sign(hB)*90)
        
    routAs, BAs = elTheory(hA, sepA)
    routBs, BBs = elTheory(hB, sepB)
    
    normAs = BAs/BA0
    normBs = BBs/BB0
    
    massA = mA/normAs
    massB = mB/normBs
    massDiff = massA - massB
    minIdx = np.where(np.abs(massDiff) == np.min(np.abs(massDiff)))[0]
    #for i in range(len(sepA)):
    #    print (i, dirs[i], massDiff[i], sepA[i], sepB[i], normAs[i], normBs[i])
    
    # Get the direction and corrected mass
    outDir = dirs[minIdx[0]]
    outMass = massA[minIdx[0]]
    return outDir, outMass

def comboFindDir(Aparams, Bparams):
    # params are [satLon, h, mass, CPA]
    lonA, hA, mA, cpaA = Aparams
    lonB, hB, mB, cpaB = Bparams
    
    # readjust h sign based on CPA
    hA = np.abs(hA)
    if cpaA < 180.:
        hA *= -1

    hB = np.abs(hB)
    if cpaB < 180.:
        hB *= -1
    
    # Standard mass calc (same as other version)
    # Get plane of sky in 0-360 deg
    posA = (lonA + 90 + 360) % 360
    posB = (lonB + 90 + 360) % 360
    if hA > 0:
        arng = [(posA-90 +360) % 360, (posA+90 +360) % 360]
    else:
        arng = [(posA-270 +360), (posA-90 + 360) % 360]
    if arng[1] < arng[0]: arng[1] += 360
        
    if hB > 0:
        brng = [(posB-90 + 360) % 360, (posB+90 + 360) % 360]
    else:
        brng = [(posB-270 +360) % 360, (posB-90 + 360) % 360]
        
    if brng[1] < brng[0]: brng[1] += 360
    
    lonMin = np.max([arng[0], brng[0]])    
    lonMax = np.min([arng[1], brng[1]])   
    if lonMax < lonMin:
        lonMax += 360 
    dirs = np.arange(lonMin,lonMax,1.)
    sepA = dirs - (lonA + np.sign(hA)*90)
    sepB = dirs - (lonB + np.sign(hB)*90)
    
    sepA[np.where(sepA > 180.)] = 360 - sepA[np.where(sepA > 180.)] 
    sepB[np.where(sepB > 180.)] = 360 - sepB[np.where(sepB > 180.)] 
    
    routA0, BA0 = elTheory(hA, 0)
    routB0, BB0 = elTheory(hB, 0)
    
    routAs, BAs = elTheory(hA, np.abs(sepA))
    routBs, BBs = elTheory(hB, np.abs(sepB))
    
    normAs = BAs/BA0
    normBs = BBs/BB0

    massA = mA/normAs
    massB = mB/normBs
    massDiff = massA - massB
    
    
    # Adding in calc of lat and deproj height
    zA = np.cos(cpaA*np.pi/180.) * np.abs(hA)
    zB = np.cos(cpaB*np.pi/180.) * np.abs(hB)
    # Non vertical component of projected radius
    projRA = np.sqrt(hA**2 - zA**2)
    projRB = np.sqrt(hB**2 - zB**2)

    # Separation between PoS and CME A
    sep1 = np.abs(lonA + 90 - dirs) % 360
    sep1[np.where(sep1 > 180.)] = 360 - sep1[np.where(sep1 > 180.)]
    sep2 = np.abs(lonA - 90 - dirs) % 360
    sep2[np.where(sep2 > 180.)] = 360 - sep2[np.where(sep2 > 180.)]
    sepA = np.array([np.min([sep1[i], sep2[i]]) for i in range(len(sep1))])
    # Separation between PoS and CME B
    sep1 = np.abs(lonB + 90 - dirs) % 360
    sep1[np.where(sep1 > 180.)] = 360 - sep1[np.where(sep1 > 180.)]
    sep2 = np.abs(lonB - 90 - dirs) % 360
    sep2[np.where(sep2 > 180.)] = 360 - sep2[np.where(sep2 > 180.)]
    sepB = np.array([np.min([sep1[i], sep2[i]]) for i in range(len(sep1))])
    # Unproject the equatorial component
    dpRA = projRA / np.cos(sepA*np.pi/180.)
    dpRB = projRB / np.cos(sepB*np.pi/180.)
    # Add back in vert components to get full R of CoM
    actRA = np.sqrt(dpRA**2 + zA**2)
    actRB = np.sqrt(dpRB**2 + zB**2)
    rDiff = actRA - actRB
    # Get latitude
    latA = np.arctan2(zA, np.abs(actRA)) * 180 / np.pi
    latB = np.arctan2(zB, np.abs(actRB)) * 180 / np.pi
    latDiff = latA - latB
    
    sclFactors = [0.5*(mA+mB), 1., 5.]
    sumScores = []
    for i in range(len(sepA)):
        sumScore = np.abs(massDiff[i]/sclFactors[0]) + np.abs(rDiff[i]/sclFactors[1]) + np.abs(latDiff[i]/sclFactors[2])
        #print (dirs[i], massDiff[i]/sclFactors[0], rDiff[i]/sclFactors[1], latDiff[i]/sclFactors[2], sumScore)
        sumScores.append(sumScore)
    bestIdx = np.where(sumScores == np.min(sumScores))[0]

    # Get the direction and corrected mass
    outDir = dirs[bestIdx[0]]
    outMassA = massA[bestIdx[0]]
    outMassB = massB[bestIdx[0]]
    outLatA = latA[bestIdx[0]]
    outLatB = latB[bestIdx[0]]
    outRA = actRA[bestIdx[0]]
    outRB = actRB[bestIdx[0]]
    
    
    return outDir, outMassA, outMassB, outLatA, outLatB, outRA, outRB


def test_it():
    # Get the satellite positions
    timeA = parse_time('2008-12-14')
    staLoc = get_horizons_coord('STEREO-A', timeA)
    stbLoc = get_horizons_coord('STEREO-B', timeA)
    lonA = staLoc.lon.degree
    lonB = stbLoc.lon.degree

    # get the difference between them
    sep  = lonA - lonB

    mA = 2.23
    mB = 2.57
    # height needs to be signed!
    rA = -10
    rB = 10

    outDir, outMass = findDir([lonA, rA, mA], [lonB, rB, mB])
    print (lonA, rA, mA, lonB, rB, mB, outDir, outMass)
    
def runCORSETcases(yr=None, sclB=1, fancyMode=False):
    # import the data structure
    f =  open('allRes.pkl', 'rb')
    res, id2time, time2time = pickle.load(f)
    f.close() 
    
    if sclB == 1:
        f2 = open('temp_CORSET_deProj.dat', 'w')
    else:
        f2 = open('CORSET_deProjSCL.dat', 'w')
        
    keys = np.sort(np.array([a for a in res.keys()]))
    counter = 0
    #for key in ['20111001_105400']:
    for key in keys:
        doIt = True
        if yr:
            if key[:4] != yr:
                doIt = False
        if doIt:
            print (key)
            myRes = res[key]
            if (myRes.CCoMtimesA[0]!=None) & (myRes.CCoMtimesB[0]!=None):
                # get matching times
                matcht = []
                idxa, idxb = [], []
                htsa, htsb = [], []
            
                for i in range(len(myRes.CCoMtimesA)):
                    aTime = myRes.CCoMtimesA[i]
                    if aTime in myRes.CCoMtimesB:
                        matcht.append(aTime)
                        idxa.append(i)
                        nowIdB = np.where(myRes.CCoMtimesB == aTime)[0][0]
                        idxb.append(nowIdB)
                        htsa.append(myRes.CCoMhtRsA[i])
                        htsb.append(myRes.CCoMhtRsB[nowIdB])
                htsa = np.array(htsa)
                htsb = np.array(htsb)
                    
                timeA = parse_time(matcht[0])
                staLoc = get_horizons_coord('STEREO-A', timeA)
                stbLoc = get_horizons_coord('STEREO-B', timeA)
                lonA = staLoc.lon.degree
                lonB = stbLoc.lon.degree
            
                # Need to pull mass and cpa (to get sign of ht)
                mAs, mBs = [], []
                cpaAs, cpaBs = [], []
                hFAs, hFBs = [], []
                
                for aTime in matcht:
                    thisIdxA = np.where(myRes.CORSETtimesA == aTime)[0]
                    mAs.append(myRes.CORSETmassesA[thisIdxA[0]])
                    cpaAs.append(myRes.CCoM_CPAsA[thisIdxA[0]])
                    hFAs.append(myRes.CCoM_hFsA[thisIdxA[0]])
                    
                    thisIdxB = np.where(myRes.CORSETtimesB == aTime)[0]
                    mBs.append(myRes.CORSETmassesB[thisIdxB[0]])
                    cpaBs.append(myRes.CCoM_CPAsB[thisIdxB[0]])
                    hFBs.append(myRes.CCoM_hFsB[thisIdxB[0]])
                cpaAs = np.array(cpaAs)
                cpaBs = np.array(cpaBs)
            
                #cpaA = myRes.CORSETcpaMA % 360
                #cpaB = myRes.CORSETcpaMB % 360
                cpaA = np.mean(cpaAs) % 360
                cpaB = np.mean(cpaBs) % 360

                # cpa defined counterclock wrt to solar N
                # h should be negative if to the left in FoV -> cpa in 0-180
                if cpaA < 180.:
                    htsa *= -1

                if cpaB < 180.:
                    htsb *= -1
                for i in range(len(idxa)):
                    #if (np.abs(htsa[i]) > 1 ) & (np.abs(htsb[i]) > 1 ):
                    #    print(key, lonA, htsa[i], mAs[i], lonB, htsb[i], mBs[i])
                    if not fancyMode:
                        try:
                            outDir, outMass = findDir([lonA, htsa[i], mAs[i]], [lonB, htsb[i], sclB*mBs[i]])
                        except:
                            outDir, outMass = None, None
                            print ('Error in ', key)
                        #print (lonA, htsa[i], mAs[i], cpaAs[i])    
                        #print (lonB, htsb[i], mBs[i], cpaBs[i])   
                        if outDir:
                            # Get the deproj lats and heights
                            # Get z component from CoM height and CPA
                            zA = np.cos(cpaAs[i]*np.pi/180.) * np.abs(htsa[i])
                            zB = np.cos(cpaBs[i]*np.pi/180.) * np.abs(htsb[i])
                            # Non vertical component of projected radius
                            projRA = np.sqrt(htsa[i]**2 - zA**2)
                            projRB = np.sqrt(htsb[i]**2 - zB**2)
                            # Separation between PoS and CME A
                            sep1 = np.abs(lonA + 90 - outDir) % 360
                            if sep1 > 180: sep1 = 360 - sep1
                            sep2 = np.abs(lonA - 90 - outDir) % 360
                            if sep2 > 180: sep2 = 360 - sep2
                            sepA = np.min([sep1, sep2])
                            # Separation between PoS and CME B
                            sep1 = np.abs(lonB + 90 - outDir) % 360
                            if sep1 > 180: sep1 = 360 - sep1
                            sep2 = np.abs(lonB - 90 - outDir) % 360
                            if sep2 > 180: sep2 = 360 - sep2
                            sepB = np.min([sep1, sep2])
                            # Unproject the equatorial component
                            dpRA = projRA / np.cos(sepA*np.pi/180.)
                            dpRB = projRB / np.cos(sepB*np.pi/180.)
                            # Add back in vert components to get full R of CoM
                            actRA = np.sqrt(dpRA**2 + zA**2)
                            actRB = np.sqrt(dpRB**2 + zB**2)
                            # Get latitude
                            latA = np.arctan2(zA, np.abs(actRA)) * 180 / np.pi
                            latB = np.arctan2(zB, np.abs(actRB)) * 180 / np.pi
                        
                            # Scale front distance same as CoM dist
                            dp_hFA = hFAs[i] * actRA / htsa[i]
                            dp_hFB = hFBs[i] * actRB / htsb[i]
                        
                            # output is key, IDA, IDB, time, lonCME, mass
                            # lonA, sepA, cpaA, latA, proj mA, projCOM htA, deprojCOM htA, p hFA, dp hFA, 
                            # lonB, sepB, cpaB, latB, proj mB, projCOM htB, deprojCOM htB, p hFB, dp hFB, 

                            outStuff = [key, myRes.CORSETidA, myRes.CORSETidB, matcht[i], '{:.2f}'.format(outDir), '{:.3f}'.format(outMass), '{:.2f}'.format(lonA), '{:.2f}'.format(sepA), '{:.2f}'.format(cpaAs[i]), '{:.2f}'.format(latA), '{:.3f}'.format(mAs[i]), '{:.3f}'.format(htsa[i]),'{:.3f}'.format(actRA),'{:.3f}'.format(hFAs[i]), '{:.3f}'.format(dp_hFA),   '{:.2f}'.format(lonB),  '{:.2f}'.format(sepB), '{:.2f}'.format(cpaBs[i]), '{:.2f}'.format(latB), '{:.3f}'.format(mBs[i]), '{:.3f}'.format(htsb[i]),'{:.3f}'.format(actRB),'{:.3f}'.format(hFBs[i]), '{:.3f}'.format(dp_hFB),]
                            outLine = ''
                            for thing in outStuff:
                                outLine += thing + ' '
                            #print (outLine)
                            f2.write(outLine+'\n')
                    else:
                        outDir, outMassA, outMassB, outLatA, outLatB, outRA, outRB = comboFindDir([lonA, htsa[i], mAs[i], cpaAs[i]], [lonB, htsb[i], sclB*mBs[i], cpaBs[i]])
                        outStuff = [key, myRes.CORSETidA, myRes.CORSETidB, matcht[i], '{:.2f}'.format(outDir), '{:.3f}'.format(outMassA), '{:.2f}'.format(lonA), '{:.2f}'.format(cpaAs[i]), '{:.2f}'.format(outLatA), '{:.3f}'.format(mAs[i]), '{:.3f}'.format(htsa[i]),'{:.3f}'.format(outRA),  '{:.3f}'.format(outMassB), '{:.2f}'.format(lonB),  '{:.2f}'.format(cpaBs[i]), '{:.2f}'.format(outLatB), '{:.3f}'.format(mBs[i]), '{:.3f}'.format(htsb[i]),'{:.3f}'.format(outRB)]
                        #print (outStuff)
                #print (sd)
    f2.close()


runCORSETcases(fancyMode=True) # sclB = 1.19

#Aparam = [6.536990972964404, -3.1349589882116558, 0.92, 94.5] #20070529_090730
#Bparam = [-3.5020973738621137, -4.013056340118309, 0.5, 84.0]

#Aparam = [67.165917585287935, 4.90940409210772, 3.02, 248.0] # 20070605_030730
#Bparam = [-3.981709828014459, 5.815735114115493, 0.57, 261.0]
#print (findDir(Aparam[:-1], Bparam[:-1]))
#print (comboFindDir(Aparam, Bparam))

#print (findDir([85.2, -10.03, 2.8], [272.7, 9.75, 1.58]))
#print (findDir([87.06, -7.29, 14.08], [266.26, 6.31, 8.06]))
#print (findDir([103.91, -6.70, 8.36], [262.24, 7.18, 6.13]))