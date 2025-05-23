import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicSpline


dtor = np.pi / 180.



def inShape(checkPt, inPt, outPt, inFit, outFit):
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
    if len(np.where(minDists < 0.1 * maxwid)[0]) / len(minDists) > 0.9:
        theInFit = lambda x: 0*x
    else:
        # Get polar generating function for inner edge
        theInFit = CubicSpline(minDangs, minDists, bc_type='natural')
    
    fakeDists = theInFit(minDangs)
    
    fig = plt.figure()
    plt.scatter(minDangs, minDists)
    plt.plot(minDangs, fakeDists,'--')
    plt.show()
    
    
    return xI, yI, theInFit
    

data = np.genfromtxt('temp.txt', dtype = float)
xs, ys = data[:,0], data[:,1]

# Get inner edge center and generating func
xI, yI, theInFit = getInner(xs, ys)

# Get outer edge center and generating func
xO, yO, theOutFit = getOuter(xs,ys)

# Make fake edges for plottin
fakeAng = np.linspace(-np.pi, np.pi,100)

fakeInDist = theInFit(fakeAng)
fakeInXs = xI + fakeInDist * np.cos(fakeAng)
fakeInYs = yI + fakeInDist * np.sin(fakeAng)

fakeDist = theOutFit(fakeAng)
fakeXs = xO + fakeDist * np.cos(fakeAng)
fakeYs = yO + fakeDist * np.sin(fakeAng)


# Loop through grid and determine if within outer bound but out of inner bound
inPt = np.array([xI, yI])
outPt = np.array([xO, yO])

nGrid = 1024
# array indexed as [y,x]
mask = np.zeros([nGrid, nGrid])-10
# get min/max range of points, no reason to check outside
minx, maxx = int(np.min(xs)), int(np.max(xs))
miny, maxy = int(np.min(ys)), int(np.max(ys))
longys = np.arange(miny,maxy+1)
for i in range(maxx - minx + 1):
    longxs = i * np.ones(len(longys)) + minx
    testPts = np.transpose([longxs,longys])
    goodys = longys[inShape(testPts,inPt, outPt, theInFit, theOutFit )]
    mask[goodys,i+minx] = 0

fig = plt.figure()
plt.contourf(mask, [-10,-3,0])
#plt.scatter(xs, ys)
plt.plot(inPt[0], inPt[1], 'ro')
plt.plot(outPt[0], outPt[1], 'ro')
plt.plot(fakeXs, fakeYs, 'go')
plt.plot(fakeInXs, fakeInYs, 'go')
plt.show()