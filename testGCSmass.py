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
    print (np.mean(dLR), critvalL, np.mean(dRL), critvalR)
    if (np.mean(dLR) < critvalL) & (np.mean(dRL) < critvalR):
        theInFit = lambda x: 0*x
        return xI, yI, theInFit
    
    fig = plt.figure()
    plt.scatter(xs[:200], ys[:200], c='r')
    plt.scatter(xs[-200:], ys[-200:], c='b')
    plt.plot(np.mean(xLR), np.mean(yLR), 'ko')
    plt.plot(np.mean(xLL), np.mean(yLL), 'ko')
    #plt.plot(xI, yI, 'ro')
    plt.show()
    
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
    
    fig = plt.figure()
    plt.plot(myAng, dist, 'bo')
    plt.plot(minDangs, minDists, 'ro')
    plt.show()
    
    # If 90% of the inner edge is less < 10% max size then assume no gap
    # These numbers are a guess at cutoffs, might need to revise
    print (len(np.where(minDists < 0.1 * maxwid)[0]) / len(minDists))
    scatterVal =  np.var(np.diff(minDists))/np.mean(minDists)
    if len(np.where(minDists < 0.1 * maxwid)[0]) / len(minDists) > 0.9:
        print ('here')
        theInFit = lambda x: 0*x
    elif scatterVal > 5:
        theInFit = lambda x: 0*x
    else:
        # Get polar generating function for inner edge
        theInFit = CubicSpline(minDangs, minDists, bc_type='natural')
        
    return xI, yI, theInFit



# Path to the corset catalog directory
#basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
basePath  = '/Users/kaycd1/STEREO_Mass/testFiles/'
# Mass files in local directory (for now?)
massPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Output path, storing with mass for now but keeping flexible
outPath = '/Users/kaycd1/STEREO_Mass/MassFits/'



#fname = 'MassFits/20120712_182400_d4c2A_mass.fts'
#lat, lon, tilt, AW, kap, ht = -15, 1.6, 60, 69, 0.36, 15

#fname = 'MassFits/20070605_083730_d4c2B_mass.fts'
#lat, lon, tilt, AW, kap, ht = -10.1,    110.5,    1.7,   27.9,  0.453,  13.57


#fname = 'MassFits/20070608_023730_d4c2B_mass.fts'
#lat, lon, tilt, AW, kap, ht = -11.7,     67.5,  -10.1,   17.0 , 0.259,  13.50

fname = 'MassFits/20100524_175400_d4c2B_mass.fts'
#lat, lon, tilt, AW, kap, ht = -2.8,     10.2,  -68.8,   31.6,  0.234,  14.29
lat, lon, tilt, AW, kap, ht = -3.9,      6.5,  -31.3,   20.1,  0.520,  11.93

#fname = 'MassFits/20071116_150754_d4c2B_mass.fts'
#lat, lon, tilt, AW, kap, ht = -14.5,    124.3 ,   5.6  , 18.4,  0.324 , 11.64


#data, header = fits.getdata(fname, header=True)
myMap = sunpy.map.Map(fname)
#data = myMap0.data
#myMap = sunpy.map.Map(data2, myMap0.meta)

# Open a mask file
#savFile = basePath +'A/2012/07/120712/2087.1938_0A/2087.1938_0A.sav'
savFile = basePath +'2087.1938_0A.sav'
sav_data = readsav(savFile)
masks = sav_data['cat'][0][10]        
        
# Get the appropriate mask from corset results
nowMask =   np.array(masks[4])

# Get the sun center and Rsun in pix
obsLon = myMap.observer_coordinate.lon.degree
obsLat = myMap.observer_coordinate.lat.degree
obsR = myMap.observer_coordinate.radius.m
skyPt = SkyCoord(obsLon*u.deg, obsLat*u.deg, 1*u.solRad,frame="heliographic_stonyhurst", obstime=myMap.date)
centS = myMap.wcs.world_to_pixel(skyPt)
sx, sy = centS[0], centS[1]

# Don't actually need rsun so much as occulted region
'''if 'rsun' in myMap.meta:
    myRs = myMap.meta['rsun'] # in arcsec
else:
    myDist = myMap.observer_coordinate.radius.m / 7e8
    myRs   = np.arctan2(1, myDist) * 206265
pixRs = myRs/myMap.scale[0].to_value() * 4'''

# took CoRe and made up some other numbers for this case
#lat, lon, tilt, AW, kap, ht = -15, 1.6, 60, 69, 0.36, 15
#lat, lon, tilt, AW, kap, ht = -8.9, 0.3, 58.1,25.4, 0.37, 15
#lat, lon, tilt, AW, kap, ht = -18, 160, 90, 45, 0.3, 15
#lat, lon, tilt, AW, kap, ht = -10.1,    110.5,    1.7,   27.9,  0.453,  13.57
pts = getGCS(lon, lat, tilt, ht, kap, AW, nleg=10)
pixXs = []
pixYs = []
for pt in pts:
    skyPt = SkyCoord(x=pt[0], y=pt[1], z=pt[2], unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
    myPt2 = myMap.world_to_pixel(skyPt)
    myX, myY = myPt2.x.to_value(), myPt2.y.to_value()
    pixXs.append(myX)
    pixYs.append(myY)
    #print (myX, myY)

pixXs = np.array(pixXs)
pixYs = np.array(pixYs)


# Make a GCS mask
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
    longxs = i * np.ones(len(longys)) + minx
    testPts = np.transpose([longxs,longys])
    goodys = longys[inShape(testPts,inPt, outPt, theInFit, theOutFit )]
    gMask[goodys,i+minx] = -3
    goodyOOs = longys[inShape(testPts,inPt, outPt, theInFit, theOutFit, outOnly=True)]
    gMaskOO[goodyOOs,i+minx] = -3
occultIdx = np.where(myMap.data == 0)
gMask[occultIdx] = 0
gMaskOO[occultIdx] = 0

print ('CORSET Mass:', np.sum(myMap.data*(nowMask == -3))/1e15)
print ('GCS Mass:   ', np.sum(myMap.data*(gMask == -3))/1e15)
print ('GCS Mass OO:', np.sum(myMap.data*(gMaskOO == -3))/1e15)

#print (gMask.shape)
#for i in range(len(pixXs)):
#    print (int(pixXs[i]),int(pixYs[i]))
#    gMask[int(pixYs[i]),int(pixXs[i])] = -3

#fig = plt.figure()
#plt.plot(pixXs, pixYs, 'o')
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection=myMap)
myMap.plot_settings['cmap'] = matplotlib.colormaps['Greys_r']
myMap.plot(axes=ax, clip_interval=(25, 99)*u.percent)
ax.contour(nowMask, levels=[-3], colors=['red'])
ax.contour(gMaskOO, levels=[-3], colors=['teal'])
ax.contour(gMask, levels=[-3], colors=['blue'])
ax.scatter(pixXs, pixYs, color='g')
ax.scatter(sx, sy, color='w')
plt.show()