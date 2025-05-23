import numpy as np
import matplotlib.pyplot as plt
import datetime 
from sunpy.coordinates import HeliocentricEarthEcliptic, get_body_heliographic_stonyhurst, get_horizons_coord
from sunpy.time import parse_time

# make label text size bigger
plt.rcParams.update({'font.size':14})

dtor = np.pi / 180.

yrs = ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']

t0 = datetime.datetime(2007,1,1)


def setupData(MassFile):
    # Build two dictionaries with everything sorted into A/B
    Adict = {}
    Bdict = {}
    data = np.genfromtxt(MassFile, dtype=str)
    for i in range(len(data[:,0])):
        myID = data[i,0]
        mySat = myID[-1]
        if mySat == 'A':
            Adict[myID] = data[i,:]
        elif mySat == 'B':
            Bdict[myID] = data[i,:]

    print ('A cases: ', len(Adict.keys()))
    print ('B cases: ', len(Bdict.keys()))
            
    # Figure out who has friends with masses            
    goodPairs = []            
    for key in Adict.keys():
        myData = Adict[key]
        myFriend = myData[-1].replace('c', '')
        if myFriend in Bdict.keys():
            goodPairs.append([key, myFriend])
    return Adict, Bdict, goodPairs

def collectThings(Adict, Bdict, goodPairs, savePairPositions = False):        
    # Collect the things for the plot
    times = []
    dts   = []
    massA = []
    massB = []
    lonAs = []
    lonBs = []
    counter = 0
    i = 0 
    if savePairPositions:
        f1 = open('pairPositions.dat', 'w')
    for pair in goodPairs:
        aVals = Adict[pair[0]]
        bVals = Bdict[pair[1]]
        i+=1
        #print (i)
        if (float(aVals[4]) > 0) & (float(bVals[4]) > 0):
            massA.append(float(aVals[4]))
            massB.append(float(bVals[4]))
            time = datetime.datetime.strptime(aVals[1], "%Y%m%dT%H%M%S" )
            dts.append(time)
            deltat = (time - t0).total_seconds()/60. / 60 / 24/ 365.
            times.append(deltat)
            # use PAs to real world coord est, not factoring in sat lat (from Mierla/Balmaceda)
            if savePairPositions:
                timeA = parse_time(aVals[1])
                staLoc = get_horizons_coord('STEREO-A', timeA)
                outstr = aVals[0] +' ' + aVals[1] +' '+ '{:6.2f}'.format((staLoc.lon.deg)) + ' ' + '{:6.2f}'.format((staLoc.lat.deg)) + ' ' + '{:6.2f}'.format((staLoc.radius.m)/7e8)
                f1.write(outstr+'\n')
                #lonAs.append(staLoc.lon.rad)
                timeB  = parse_time(bVals[1])
                stbLoc = get_horizons_coord('STEREO-B', timeB)
                outstr = bVals[0]+' ' + bVals[1]+ ' ' + '{:6.2f}'.format((stbLoc.lon.deg)) + ' ' + '{:6.2f}'.format((stbLoc.lat.deg)) + ' ' + '{:6.2f}'.format((stbLoc.radius.m)/7e8)
                #lonBs.append(stbLoc.lon.rad) 
                f1.write(outstr+'\n')            
        else: 
            counter += 1
    if savePairPositions:
        f1.close() 
        return times, dts, massA, massB, counter, lonsA, lonsBs
    else:      
        return times, dts, massA, massB, counter

def calcStats(times, dts, massA, massB, counter):        
    massA = np.array(massA)
    massB = np.array(massB)
    times = np.array(times)
    dts   = np.array(dts)

    avgMass  = (massA+massB) * 0.5
    massDiff = (massA - massB)/ avgMass
    sclA = massA / avgMass
    sclB = massB / avgMass
    massDiffA = (massA - avgMass)/ avgMass
    massDiffB = (massB - avgMass)/ avgMass

    # Data analysis
    print (len(dts), 'matched cases')
    print (counter, 'cases with negative masses')
    print ('')
    print ('           Mean  Med')
    print ('A/B:      ', '{:.2f}'.format(np.mean(massA / massB)), '{:.2f}'.format(np.median(massA / massB)))
    print ('A/avg:    ', '{:.2f}'.format(np.mean(sclA)), '{:.2f}'.format(np.median(sclA)))
    print ('B/avg:    ', '{:.2f}'.format(np.mean(sclB)), '{:.2f}'.format(np.median(sclB)))
    print ('|A-B|/avg:', '{:.2f}'.format(np.mean(np.abs(massDiff))), '{:.2f}'.format(np.median(np.abs(massDiff))))
    print (len(massDiff), len(np.where(np.abs(massDiff)<1)[0]), len(np.where(np.abs(massDiff)<0.5)[0]))

# Scatter A vs B
if False:
    fig, axes = plt.subplots(1, 1, figsize=(8,8))
    axes.plot([0,25], [0,25], 'k--', lw=3)
    c = axes.scatter(massA, massB, c=times, cmap='plasma')
    axes.set_aspect('equal')
    ax0pos = axes.get_position()
    fig.subplots_adjust(top=0.9)
    cbar_ax = fig.add_axes([ax0pos.x0, 0.94, ax0pos.width, 0.02])
    cbar = fig.colorbar(c, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_title('Years since 2007')   
    axes.set_xlabel('STEREO A Mass (10$^{15}$ g)')     
    axes.set_ylabel('STEREO B Mass (10$^{15}$ g)')     
    #axes.set_xscale('log')
    #axes.set_yscale('log')
    plt.savefig('fig_maxMassCompLog.png')
    
if False:
    fig, axes = plt.subplots(1, 1, figsize=(8,6.5))
    axes.plot([t0,t0+datetime.timedelta(days=(8*365))], [0,0], 'k-', zorder=0)
    axes.plot([t0,t0+datetime.timedelta(days=(8*365))], [1,1], 'k:', zorder=0)
    axes.plot([t0,t0+datetime.timedelta(days=(8*365))], [-1,-1], 'k:', zorder=0)
    c = plt.scatter(dts, massDiff, c=massA, cmap='plasma',vmin=0, vmax=15)
    ax0pos = axes.get_position()
    fig.subplots_adjust(top=0.85)
    cbar_ax = fig.add_axes([ax0pos.x0, 0.9, ax0pos.width, 0.02])
    cbar = fig.colorbar(c, cax=cbar_ax, orientation='horizontal')
    axes.set_ylim(-2,2)
    axes.set_xlim(t0,t0+datetime.timedelta(days=(8*365)))
    cbar.ax.set_title('Avg. CME Mass (10$^{15}$ g)')   
    axes.set_ylabel('(Mass$_A$ - Mass$_B$) / Avg. Mass')     
    axes.set_xlabel('Time')
    plt.savefig('fig_maxMassTimeline.png')
    
if False:
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10,5))
    bins = np.linspace(-1,1,20)
    axes[0].hist((massA-avgMass)/avgMass, bins=bins)
    axes[1].hist((massB-avgMass)/avgMass, bins=bins)
    ylim = axes[0].get_ylim()
    axes[0].plot([0,0],ylim, 'k--')
    axes[1].plot([0,0],ylim, 'k--')
    axes[0].set_ylim(ylim)
    axes[0].set_ylabel('Counts')
    axes[0].set_xlabel('(Mass$_A$ - Avg. Mass) / Avg. Mass')
    axes[1].set_xlabel('(Mass$_B$ - Avg. Mass) / Avg. Mass')
    plt.savefig('fig_maxMassHistos.png')
    
if False:
    CMElons = []
    CMElats = []
    CMErs   = []
    # we are assuming these obs are close enough to same time in A/B, if heights are diff
    # this will all be thrown off
    for i in range(len(lonAs)):
        pair = goodPairs[i]
        aVals = Adict[pair[0]]
        bVals = Bdict[pair[1]]
        lonA = lonAs[i]
        lonB = lonBs[i]
        cpaA = float(aVals[5])
        rA   = float(aVals[3])
        cpaB = float(bVals[5])
        rB   = float(bVals[3])
        
        x = (rA * np.sin(cpaA) * np.cos(lonB) - rB * np.sin(cpaB) * np.cos(lonA)) 
        # 99% certain the - in Laura's eq 4 should be in numerator first term, not on whole thing
        # the calc matches the Mierla version if we correct as such
        y = -(rA * np.sin(cpaA) * np.sin(lonB) - rB * np.sin(cpaB) * np.sin(lonA)) 
        lon = np.arctan(y/x) / dtor
        z = 0.5*(rB * np.cos(cpaB) + rA * np.cos(cpaA))
        lat = np.arctan(z/np.sqrt(x**2 + y**2))/dtor
        r3d = np.sqrt(x**2 + y**2 + z**2)
        
        CMElons.append(lon)
        CMElats.append(lat)
        CMErs.append(r3d)
    
    CMElons = np.array(CMElons)
    fig, axes = plt.subplots(1, 1, figsize=(8,6.5))
    axes.scatter(CMElons, massDiff)
    plt.show()
    
Adict, Bdict, goodPairs = setupData('CORSET_MaxMasses.dat')
times, dts, massA, massB, counter = collectThings(Adict, Bdict, goodPairs, savePairPositions=False)
calcStats(times, dts, massA, massB, counter)