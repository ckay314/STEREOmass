import numpy as np
import sunpy.map
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from scipy.io import readsav
import sys, os
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicSpline

sys.path.append(os.path.abspath('/Users/kaycd1/wombat/')) 
from pyGCS import getGCS

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)


global basePath, massPath, outPath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
# Path to the fits files
fitsPathA = '/Volumes/SRP/vourla1/secchi/lz/L0/a/img/cor2/'
fitsPathB = '/Volumes/SRP/vourla1/secchi/lz/L0/b/img/cor2/'
# Mass files in local directory (for now?)
massPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Output path, storing with mass for now but keeping flexible
outPath = '/Users/kaycd1/STEREO_Mass/GCSimages/'



# Helper functions 
def inShape(checkPt, inPt, outPt, inFit, outFit, outOnly = False):
    if len(checkPt.shape) == 1:
       isIn = np.array([True])
       checkPt = np.array([checkPt])
    elif len(checkPt.shape) == 2:
        isIn = np.array([True for i in range(len(checkPt))])
    
    # get polar coords wrt in/outPts
    dxI  = checkPt[:,0] - inPt[0]
    dyI  = checkPt[:,1] - inPt[1]
    rI   = np.sqrt(dxI**2 + dyI**2)
    angI = np.arctan2(dyI, dxI)

    dxO = checkPt[:,0] - outPt[0]
    dyO = checkPt[:,1] - outPt[1]
    rO  = np.sqrt(dxO**2 + dyO**2)
    angO = np.arctan2(dyO, dxO)
    
    # check if outside of inner bound
    if not outOnly:
        critRI = inFit(angI)
        inIdx = np.where(rI < critRI)[0]
        isIn[inIdx] = False
    
    # check if inside of outer bounds    
    critRO = outFit(angO)
    outIdx = np.where(rO > critRO)[0]
    isIn[outIdx] = False
    
    return isIn

def getOuter(xs, ys):
    # Get outer edge of set xy
    hull = ConvexHull(np.transpose(np.array([xs,ys])))
    outs = hull.vertices
    
    # Select the edge points
    exs, eys = xs[outs], ys[outs]
    
    # Get the mean of the edges
    xO, yO = np.mean(exs), np.mean(eys)
    
    # Need to check if there are gaps where it wasn't convex
    dists = np.zeros(len(exs))
    dists[1:-1] = np.sqrt((exs[2:] -exs[:-2])**2 +  (eys[2:] - eys[:-2])**2)
    dists[0] = np.sqrt((exs[1] -exs[-1])**2 + (eys[1] - eys[-1])**2)
    dists[-1] = np.sqrt((exs[0] -exs[-2])**2 + (eys[0] - eys[-2])**2)
    meanDist, stdDist = np.mean(dists), np.std(dists)
    
    # Pull the points with large neighbor dists
    gapEdges = np.where(dists > (meanDist + 2 * stdDist))[0]
    
    # Check if adjacent pairs
    gapPairs = []
    for idx in gapEdges:
        if idx + 1 in gapEdges:
            gapPairs.append([idx, idx+1])
        elif idx == 0:
            if len(exs) - 1 in gapEdges:
                gapPairs.append([len(exs)-1, 0])

    for aPair in gapPairs:
        idx1, idx2 = aPair
        ex1, ey1 = exs[idx1], eys[idx1]
        ex2, ey2 = exs[idx2], eys[idx2]
        thisDist = np.sqrt((ex1-ex2)**2 + (ey1-ey2)**2)
        npts = int(thisDist/meanDist*2)
        dx = (ex2 - ex1) / (npts+1)
        dy = (ey2 - ey1) / (npts+1)
        newxs = ex1 + dx * (1 + np.arange(npts))
        newys = ey1 + dy * (1 + np.arange(npts))
        exs = np.append(exs, newxs)
        eys = np.append(eys, newys)
    
    # Calc dists from midpoint
    distx, disty = exs-xO, eys-yO
    dist = np.sqrt(distx**2 + disty**2)
    myAng = np.arctan2(disty, distx) 
    angSort = np.argsort(myAng)
    dist = dist[angSort]
    myAng = myAng[angSort]

    # Make polar function for dist to edge from center as function of angle
    theOutFit = CubicSpline(myAng, dist, bc_type='natural')
    
    return xO, yO, theOutFit
    
def getInner(xs, ys):
    # Get mean of all points
    xI, yI = np.mean(xs), np.mean(ys)
    xM, yM = 0.5*(np.max(xs) + np.min(xs)), 0.5*(np.max(ys) + np.min(ys))
    
    # Extract the leg portions (this is assuming nLeg set to 10 and rest at default)
    xLR, yLR = xs[:200], ys[:200]
    xLL, yLL = xs[-200:], ys[-200:]
    
    # Get avg dist from own center and other leg center
    xLRc, yLRc = np.mean(xs[:200]), np.mean(ys[:200])
    xLLc, yLLc = np.mean(xs[-200:]), np.mean(ys[-200:])
    dLL = np.sqrt((xLR- xLRc)**2 + (yLR - yLRc)**2)
    dLR = np.sqrt((xLR- xLLc)**2 + (yLR - yLLc)**2)
    dRL = np.sqrt((xLL- xLRc)**2 + (yLL - yLRc)**2)
    dRR = np.sqrt((xLL- xLLc)**2 + (yLL - yLLc)**2)
    critvalL = 0.5*(np.mean(dLL) + np.max(dLL))
    critvalR = 0.5*(np.mean(dRR) + np.max(dRR))
    
    # If dists between less than dists from other then return
    # as no gap case
    if (np.mean(dLR) < critvalL) & (np.mean(dRL) < critvalR):
        theInFit = lambda x: 0*x
        return xI, yI, theInFit
    
    # Get width of points to check gap against
    xwid = np.max(xs) - np.min(xs)
    ywid = np.max(ys) - np.min(ys)
    maxwid = np.max([xwid,ywid])
    
    # Calc polar coords
    distx, disty = xs-xI, ys-yI
    dist = np.sqrt(distx**2 + disty**2)
    myAng = np.arctan2(disty, distx) 
    
    # Sort by angle so spline is happy
    angSort = np.argsort(myAng)
    dist = dist[angSort]
    myAng = myAng[angSort]
        
    # Make angle bins
    fakeAng = np.linspace(-np.pi, np.pi,30)
    dAng = 0.5 * (fakeAng[1] - fakeAng[0])
    minDists = []
    minDangs = []
    # Loop through the angle bins and pick the minimum r for that bin
    for i in range(len(fakeAng)):
        theseAng = np.where(np.abs(myAng - fakeAng[i]) < dAng)[0]
        if len(theseAng) > 0:
            minDists.append(np.min(dist[theseAng]))
            minDangs.append(fakeAng[i])
    minDists = np.array(minDists)
    minDangs = np.array(minDangs)      
    
    
    # If 90% of the inner edge is less < 10% max size then assume no gap
    # These numbers are a guess at cutoffs, might need to revise
    scatterVal =  np.var(np.diff(minDists))/np.mean(minDists)
    if len(np.where(minDists < 0.1 * maxwid)[0]) / len(minDists) > 0.9:
        theInFit = lambda x: 0*x
    elif scatterVal > 5:
        theInFit = lambda x: 0*x
    else:
        # Get polar generating function for inner edge
        theInFit = CubicSpline(minDangs, minDists, bc_type='natural')        
    return xI, yI, theInFit

    
def makeCORSETdict():
    data = np.genfromtxt('corsetDirectory.txt', dtype = str)
    CORdict = {}
    for i in range(len(data[:,0])):
        CORdict[data[i,0]] = data[i,1]
    return CORdict

def getCORSET(CORid, CORpath):
    savFile = CORpath+CORid+'/'+CORid+'.sav'
    infoFile = CORpath+CORid+'/'+CORid+'.info'
    # make sure it exists
    if os.path.exists(savFile):   
        # Get the masks
        try:
            sav_data = readsav(savFile)
        except:
            return None, None
        masks = sav_data['cat'][0][10]
        
        # Get the time stamps of the masks
        f1 = open(infoFile, 'r')
        infos = []
        baseIdx  = -1
        fileIdx0 = -1
        counter = -1
        # Figure out the path string we want to rm
        if 'A' in CORid:
            rmPath = fitsPathA
        elif 'B' in CORid:
            rmPath = fitsPathB
        else:
            sys.exit('Error in processing CORid')
        
        for x in f1:
            counter += 1
            if 'Base image' in x:
                baseIdx = counter +1
            elif 'Files used' in x:
                fileIdx0 = counter + 1
            infos.append(x.replace('\n', '').replace(rmPath, ''))
        f1.close
        CORfiles = infos[fileIdx0:]
        return masks, np.array(CORfiles)
    else:
        return None, None
        
def makeGCSmask(lat, lon, tilt, AW, kap, ht, myMap):
    pts = getGCS(lon, lat, tilt, ht, kap, AW, nleg=10)
    pixXs = []
    pixYs = []
    for pt in pts:
        skyPt = SkyCoord(x=pt[0], y=pt[1], z=pt[2], unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
        myPt2 = myMap.world_to_pixel(skyPt)
        myX, myY = myPt2.x.to_value(), myPt2.y.to_value()
        pixXs.append(myX)
        pixYs.append(myY)

    pixXs = np.array(pixXs)
    pixYs = np.array(pixYs)
    
    
    # Make a GCS mask
    nx = myMap.data.shape[1]
    gMask = np.zeros(myMap.data.shape)
    gMaskOO = np.zeros(myMap.data.shape)

    # Get inner edge center and generating func
    xI, yI, theInFit = getInner(pixXs, pixYs)

    # Get outer edge center and generating func
    xO, yO, theOutFit = getOuter(pixXs, pixYs)

    # Loop through and fill the GCS mask
    # get min/max range of points, no reason to check outside
    minx, maxx = int(np.min(pixXs)), int(np.max(pixXs))
    miny, maxy = int(np.min(pixYs)), int(np.max(pixYs))
    inPt = np.array([xI, yI])
    outPt = np.array([xO, yO])
    longys = np.arange(miny,maxy+1)
    for i in range(maxx - minx + 1):
        if (i+minx < nx-1) & (i+minx >= 0):
            longxs = i * np.ones(len(longys)) + minx
            testPts = np.transpose([longxs,longys])
            goodys = longys[inShape(testPts,inPt, outPt, theInFit, theOutFit )]
            gMask[goodys,i+minx] = -3
            goodyOOs = longys[inShape(testPts,inPt, outPt, theInFit, theOutFit, outOnly=True)]
            gMaskOO[goodyOOs,i+minx] = -3
    occultIdx = np.where(myMap.data == 0)
    gMask[occultIdx] = 0
    gMaskOO[occultIdx] = 0    
    return gMask, gMaskOO
    
def makeOutputs(myMap, myCORmask, gMask, gMaskOO, imName, im0Name, usrName, CORid=None, doFig=True): 
    if np.any(myCORmask):
        CORmass = '{:4.2f}'.format(np.sum(myMap.data*(myCORmask == -3))/1e15)
        nCpix = np.sum(myCORmask == -3)
    else:
        CORmass = 'None'
        nCpix   = 'None' 
    
    doGCS = False
    if type(gMask) != type(None):
        doGCS = True
        GCSmass = '{:4.2f}'.format(np.sum(myMap.data*(gMask == -3))/1e15)
        nGpix = np.sum(gMask == -3)
        GCSmassOO = '{:4.2f}'.format(np.sum(myMap.data*(gMaskOO == -3))/1e15)
        nGpixOO = np.sum(gMaskOO == -3)
        #print ('CORSET Mass:', CORmass, nCpix)
        #print ('GCS Mass:   ', GCSmass, nGpix)
        #print ('GCS Mass OO:', GCSmassOO, nGpixOO)
                
    if doFig:            
        fig = plt.figure(figsize = (8,8), frameon=False)
        ax = fig.add_subplot(projection=myMap)
        myMap.plot_settings['cmap'] = matplotlib.colormaps['Greys_r']
    
        # try and mask the occulter
        occMask = myMap.data * 0.
        occIds = np.where(myMap.data == 0)
        occMask[occIds] = -3
        
        myMap.plot(axes=ax, clip_interval=(25, 99)*u.percent)
        if np.any(myCORmask):
            ax.contour(myCORmask, levels=[-3], colors=['red'], linewidths=3)
        if doGCS:
            ax.contour(gMaskOO, levels=[-3], colors=['teal'])
            ax.contour(gMask, levels=[-3], colors=['blue'], linewidths=3)
        ax.contourf(occMask, levels=[-5,-3], colors=['black'])
        ax.set_axis_off()
        ax.set_title('')
    
        # Add im/base time and total mass as text
        ax.text(0.01, 0.06, 'Mass (10$^{15}$ g)', transform=ax.transAxes, color='w', horizontalalignment='left', verticalalignment='bottom', fontsize=16)
        if doGCS:
            if GCSmassOO != GCSmass:
                ax.text(0.01, 0.03, 'GCS: ' + GCSmass+' ('+ GCSmassOO +')', transform=ax.transAxes, color='w', horizontalalignment='left', verticalalignment='bottom', fontsize=16)
            else:
                ax.text(0.01, 0.03, 'GCS: ' + GCSmass, transform=ax.transAxes, color='w', horizontalalignment='left', verticalalignment='bottom', fontsize=16)
        if np.any(myCORmask):
            ax.text(0.01, 0.0, 'CORSET: '+CORmass, transform=ax.transAxes, color='w', horizontalalignment='left', verticalalignment='bottom', fontsize=16)
        subImName = imName.replace('.CK', '').replace('_d4c2A', '').replace('_d4c2B', '')
        subImName0 = im0Name.replace('.CK', '').replace('_d4c2A', '').replace('_d4c2B', '')
        ax.text(0.99, 0.03, subImName, transform=ax.transAxes, color='w', horizontalalignment='right', verticalalignment='bottom', fontsize=16)
        ax.text(0.99, 0.0, subImName0, transform=ax.transAxes, color='w', horizontalalignment='right', verticalalignment='bottom', fontsize=16)
        plt.subplots_adjust(hspace=0.1,left=0,right=1,top=1,bottom=0)
        
        outFile = outPath+imName
        if CORid not in [None, 'None']:
            outFile = outFile+'_'+CORid
        outFile = outFile + '_' + usrName + '_mass.png'   
        plt.savefig(outFile)
        plt.close()
    
    if doGCS:
        outstr = imName.ljust(22) + CORid.rjust(13) + CORmass.rjust(6) + str(nCpix).rjust(10) + usrName.rjust(15) + GCSmass.rjust(6) + str(nGpix).rjust(10) + GCSmassOO.rjust(6) + str(nGpixOO).rjust(10)
        print (outstr)
   
        return outstr

def processCase(aIm, aIm0, aUsr, myGCS, CORmasks, CORfiles, CORid):
    # Figure out if we have a matching CORSET case
    myCORmask = None
    if np.any(CORmasks):
        if aIm in CORfiles:
            myCORidx = np.where(CORfiles == aIm)[0]
            myCORmask = np.array(CORmasks[myCORidx[-1]])
        else:
            # Should have it but not finding for some reason
            print ('Error in finding matching CORSET mask')
    
    # Get the appropriate mass file
    slashIdx = np.char.find(aIm, '/')
    shortName = aIm[slashIdx+1:].replace('.fts', '_massPOS.fts')
    shortestName = aIm[slashIdx+1:].replace('.fts', '')
    slashIdx0 = np.char.find(aIm0, '/')
    shortName0 = aIm0[slashIdx0+1:].replace('.fts', '_massPOS.fts')
    shortestName0 = aIm0[slashIdx0+1:].replace('.fts', '')
    
    if os.path.exists(massPath+shortName.replace('.CK','')):
        massFts = massPath+shortName
    else:
        print('Cannot find mass file: ' + massPath+shortName)
        return None
    
    # Set up map for the background image 
    myMap = sunpy.map.Map(massFts)
    
    
    
    # Use map and GCS to make masks
    # OO = outer only, ignores inner leg gap in GCS 
    lat, lon, tilt, AW, kap, ht = myGCS
    gMask, gMaskOO = makeGCSmask(lat, lon, tilt, AW, kap, ht, myMap)
    
    # Make all the outputs (numbers and fig)
    outstr = makeOutputs(myMap, myCORmask, gMask, gMaskOO, shortestName, shortestName0, aUsr, CORid=CORid)
    return outstr

def runAllCases(saveIt=False):
    # Open the files with info about GCSs
    dataCore = np.genfromtxt('coreGCScases.txt', dtype=str)
    dataIms = np.genfromtxt('haveMassImage.txt', dtype=str)
    
    if saveIt:
        f1 = open('GCSmassComparisonPOS.txt', 'w')

    # Get list ocf unique COR time stamps
    timeIDs = np.unique(dataCore[:,0])
    counter = 0
    n2do    = len(timeIDs)
    
    # use to find indiv cases
    #for i in range(len(timeIDs)):
    #    print(i, timeIDs[i])
    #print(sd)
    #while counter < n2do:
    for counter in [221]:
        counter += 1
        aTime = timeIDs[counter-1]
        print ('On event ', counter, ' out of ', n2do, aTime)
                
        myIdx  = np.where(dataCore[:,0] == aTime)[0] # same as np.where(dataIms[:,0] == aTime)[0])
        myUsrs = dataIms[myIdx, 1]
        # Take out Majumdar cases bc COR1 (for now at least)
        if 'Majumdar' in myUsrs:
            rmIdx  = np.where(myUsrs == 'Majumdar')[0]
            myIdx  = np.delete(myIdx, rmIdx)
            myUsrs = np.delete(myUsrs, rmIdx)
            
        # Take out Temmer21 cases bc no height (for now at least)
        if 'Temmer21' in myUsrs:
            rmIdx  = np.where(myUsrs == 'Temmer21')[0]
            myIdx  = np.delete(myIdx, rmIdx)
            myUsrs = np.delete(myUsrs, rmIdx)

        myAims  = dataIms[myIdx, 7]
        myBims  = dataIms[myIdx, 12]
        myA0ims = dataIms[myIdx, 6]
        myB0ims = dataIms[myIdx, 11]

        # make sure we have an A event image    
        if 'None' in myAims:
            rmIdxA  = np.where(myAims == 'None')[0]
            myIdxA  = np.delete(myIdx, rmIdxA)
            myUsrsA = np.delete(myUsrs, rmIdxA)
            myAims  = np.delete(myAims, rmIdxA)
            myA0ims = np.delete(myA0ims, rmIdxA)
        else:
            myIdxA, myUsrsA = myIdx, myUsrs

        # make sure we have an A event image    
        if 'None' in myBims:
            rmIdxB  = np.where(myBims == 'None')[0]
            myIdxB  = np.delete(myIdx, rmIdxB)
            myUsrsB = np.delete(myUsrs, rmIdxB)
            myBims  = np.delete(myBims, rmIdxB)
            myB0ims = np.delete(myB0ims, rmIdxB)
        else:
            myIdxB, myUsrsB = myIdx, myUsrs

    
        if np.any(myAims):            
            # Get the corset ID (should be same for all)
            CORid = dataIms[myIdxA[0], 4]
            CORmasks = None
            CORfiles = None
            if CORid != 'None':
                # Pull the file path from the dictionary
                CORpath  = CORdict[CORid]
                CORmasks, CORfiles = getCORSET(CORid, CORpath)
                
            for i in range(len(myAims)):
                aIm = myAims[i].replace('.CK', '.fts')
                aIm0 = myA0ims[i].replace('.CK', '.fts')
                aUsr = myUsrsA[i]
                bigIdx = myIdxA[i]
                
                # Pull the GCS fit for this case
                myLine = dataCore[bigIdx, :]
                lat, lon, tilt, AW, kap, ht = float(myLine[4]), float(myLine[5]), float(myLine[6]), float(myLine[7]), float(myLine[8]), float(myLine[9])
                myGCS = [lat, lon, tilt, AW, kap, ht]
                
                outstr = processCase(aIm, aIm0, aUsr, myGCS, CORmasks, CORfiles, CORid)
                if outstr and saveIt:
                    f1.write(outstr+'\n')

        if np.any(myBims):            
            # Get the corset ID (should be same for all)
            CORid = dataIms[myIdxB[0], 9]
            CORmasks = None
            if CORid != 'None':
                # Pull the file path from the dictionary
                CORpath  = CORdict[CORid]
                CORmasks, CORfiles = getCORSET(CORid, CORpath)
                
            for i in range(len(myBims)):
                bIm = myBims[i].replace('.CK', '.fts')
                bIm0 = myB0ims[i].replace('.CK', '.fts')
                bUsr = myUsrsB[i]
                bigIdx = myIdxB[i]
                
                # Pull the GCS fit for this case
                myLine = dataCore[bigIdx, :]
                lat, lon, tilt, AW, kap, ht = float(myLine[4]), float(myLine[5]), float(myLine[6]), float(myLine[7]), float(myLine[8]), float(myLine[9])
                myGCS = [lat, lon, tilt, AW, kap, ht]
                
                outstr = processCase(bIm, bIm0, bUsr, myGCS, CORmasks, CORfiles, CORid)
                if outstr and saveIt:
                    f1.write(outstr+'\n')
    if saveIt:
        f1.close()

def justCORSET(CORid, time, time0):
    # Pull the file path from the dictionary
    CORpath  = CORdict[CORid]
    CORmasks, CORfiles = getCORSET(CORid, CORpath)
    # split up CORfiles
    if 'A' in CORid:
        subNames = np.array([item.split('/')[1].replace('_d4c2A.fts','') for item in CORfiles])
    else:
        subNames = np.array([item.split('/')[1].replace('_d4c2B.fts','') for item in CORfiles])
    if time not in subNames:
        print('Cannot find ', time, ' in CORSET masks. Options are:')
        print(subNames)
    else:   
        if 'A' in CORid: 
            ftsFile = CORpath+CORid+'/fts/'+CORid+'_'+time+'_d4c2A_mass.fts'
        else:
            ftsFile = CORpath+CORid+'/fts/'+CORid+'_'+time+'_d4c2B_mass.fts'
        # Set up map for the background image 
        myMap = sunpy.map.Map(ftsFile)
        
        idx = np.where(subNames == time)[0]
        myCORmask = CORmasks[idx[0]]
        makeOutputs(myMap, myCORmask, None, None, time, time0, 'COR', CORid=CORid)
        
    
        
    
    
    

if __name__ == '__main__':
    global CORdict
    CORdict = makeCORSETdict()
    #runAllCases(saveIt=False)
    
    # Bad example cases for the paper appendix
    CORid = '1663.5167_0A'
    time  = '20110516_023900'
    time  = '20110516_045400'
    time0 = '20110516_003900'
    
    CORid = '1663.5583_0B'
    time  = '20110516_023900'
    time  = '20110516_045400'
    time0 = '20110516_012422'
    
    
    justCORSET(CORid, time, time0)